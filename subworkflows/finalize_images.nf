#!/usr/bin/env nextflow

include { finalize; rescale; index_cellcrops; index_cellcrops_csv; cellcrops } from '../processes/finalize.nf'
include { measure_intensity; stage_as_plate } from '../processes/intensity.nf'
include { index_images as index_unscaled; index_images as index_scaled } from '../processes/staging.nf'
include { estimate_scaling_factors } from './scaling_factors.nf'
include { asBoolean } from '../lib/utils.nf'

// A path-type process output emitted via a glob (e.g. path("*.tiff")) comes through as a
// bare Path when exactly one file matches, and only as a List when more than one file
// matches (same quirk already worked around elsewhere - see the instanceof List checks in
// flatfield_estimation.nf and finalize.nf). Normalize to "the well's own output directory"
// regardless of match count.
def wellOutputDir(pathOrList) {
    return (pathOrList instanceof List) ? pathOrList[0].parent : pathOrList.parent
}

// Build a deterministic (plate, fingerprint) pair from a finalize/rescale processed_output
// channel. Nextflow's own cache hash for a directory `path` input only looks at that
// directory's own metadata, not its contents (confirmed empirically), so relying on the
// published directory path alone to detect real image changes silently misses changes
// made several directories deep.
//
// One entry per WELL, not per file: every channel file for a well is written together,
// atomically, into the same well.relpath leaf directory by a single finalize/rescale task,
// so that directory's own mtime is a sufficient stand-in for "did this well's images
// change" - confirmed empirically that Nextflow's publishDir (both explicit mode:"copy"
// and the default symlink mode this pipeline actually uses) bumps the containing
// directory's mtime when it replaces an existing published file, so an in-place
// re-publish (e.g. a well re-rescaled with new scaling factors) is still detected.
// This keeps the cost to one stat() per well regardless of channel count, instead of one
// per file.
def buildPlateFingerprint(processed_output) {
    return processed_output
        .map{ row -> tuple(row[0].plate, "${row[0].key}:${wellOutputDir(row[1]).lastModified()}") }
        .groupTuple(by: 0)
        .map{ row -> tuple(row[0], row[1].sort()) }
        .multiMap{ row -> plate: row[0]; fingerprint: row[1] }
}

// Same idea as buildPlateFingerprint, but index_cellcrops aggregates every well across
// every plate into a single output (no per-plate parameter to key off), so there's just
// one fingerprint for the whole run instead of one per plate. .collect() plays the same
// role groupTuple did above: it doesn't emit until the source channel fully closes (same
// "wait for every well" barrier timing as the .last() it replaces), but the payload is a
// deterministic, sorted list of well/mtime pairs instead of an arbitrary last-completed item.
def buildCellcropsFingerprint(h5_out) {
    return h5_out
        .map{ row -> "${row[0].key}:${wellOutputDir(row[1]).lastModified()}" }
        .collect()
        .map{ entries -> entries.sort() }
}

// Finalizes and caches images for feature extraction (rn_cache_images).
//
// Normal flow (sc_scale_in_finalize=false): finalize writes unscaled images to
// rr__processed_images/unscaled, intensities are always measured on those (even if no scaling is
// requested, to support a future QC step), then - if scaling is requested - factors are
// estimated from those same unscaled images and a separate rescale pass produces
// rr__processed_images/scaled.
//
// Fast path (sc_scale_in_finalize=true, requires sc_manualscale - enforced in checks.nf):
// scaling factors are already known upfront, so finalize applies them directly and writes
// straight to rr__processed_images/scaled, skipping the unscaled artifact, measurement, and the
// separate rescale pass entirely.
workflow finalize_images {

    take:
        finalize_in
        image_dir_file
        flatfield_out
        manifest_registration
        sc_scale_in_finalize
        sc_autoscale
        sc_manualscale
        sc_scale_slope
        sc_scale_bias
        sc_publish_unscaled
        rn_make_cellcrops
        rn_manifest_registration
        rn_publish_dir
        // Ingredients for estimating autoscale factors from the measured unscaled images
        // (only used when sc_autoscale is true)
        sc_control_list
        sc_channel_map

    main:

        // Determine if scaling needs to be run at all
        run_scaling = sc_manualscale != null || sc_autoscale

        //------------------------------------------------------------------------
        // Block 1: finalize unscaled OR direct manual scaled images wihtout caching them in between
        //------------------------------------------------------------------------
        if (sc_scale_in_finalize && run_scaling) {
            finalize_scaling_in = Channel.value(file(sc_manualscale))
            finalize_slope_in   = (sc_scale_slope != null) ? Channel.value(file(sc_scale_slope)) : Channel.value(file("NO_SLOPE"))
            finalize_bias_in    = (sc_scale_bias != null)  ? Channel.value(file(sc_scale_bias))  : Channel.value(file("NO_BIAS"))
        } else {
            finalize_scaling_in = Channel.value(file("NO_SCALE"))
            finalize_slope_in   = Channel.value(file("NO_SLOPE"))
            finalize_bias_in    = Channel.value(file("NO_BIAS"))
        }

        finalize_out = finalize(finalize_in,
                    image_dir_file,
                    flatfield_out,
                    finalize_scaling_in,
                    finalize_slope_in,
                    finalize_bias_in)

        finalize_subdir = sc_scale_in_finalize ? "scaled" : "unscaled"

        // Index whatever finalize actually produced (index_unscaled despite the name -
        // when sc_scale_in_finalize writes directly to `scaled`, this call indexes that
        // instead; it's just the alias used for "whatever finalize wrote" since a process
        // can only be invoked once per workflow scope)
        //
        // previous_completed carries a per-plate content fingerprint (see
        // buildPlateFingerprint above) instead of a throwaway barrier value - it blocks
        // until every well has finalized, same as before, but the actual value now
        // reflects whether that plate's images genuinely changed, so this only reruns
        // when it needs to instead of on nearly every run (or, worse, never noticing a
        // real change - see buildPlateFingerprint's comment for why).
        finalize_fp = buildPlateFingerprint(finalize_out.processed_output)
        if (finalize_subdir != "unscaled" || asBoolean(sc_publish_unscaled)) {
            index_unscaled(finalize_fp.fingerprint,
                            "rr__processed_images/${finalize_subdir}",
                            file(rn_publish_dir + "/rr__processed_images/${finalize_subdir}"),
                            finalize_fp.plate)
        }

        // Always measure these images, even when no scaling is requested, to
        // support a future QC step
        measure_intensity(finalize_out.processed_output, finalize_subdir)

        //------------------------------------------------------------------------
        // Block 2: default path, use cached images to determine scaling factors
        //------------------------------------------------------------------------
        if (!sc_scale_in_finalize && run_scaling) {
            if (sc_autoscale) {
                // Auto scaling: aggregate per-well parquet output into one directory per
                // plate before handing off to calculate_scaling_factors, so Nextflow
                // doesn't have to stage a well-level filelist that can run into the
                // thousands of files for a large plate
                plate_measurements = stage_as_plate(
                    measure_intensity.out.measurements
                        .map{row -> tuple(row[0].plate, row[1])}
                        .groupTuple(by: 0)
                        // groupTuple collects each well's [image, object] parquet pair into
                        // its own nested list - flatten so stage_as_plate's path input sees
                        // one flat filelist per plate, not a list-of-lists
                        .map{row -> tuple(row[0], row[1].flatten())}
                )

                // estimate_scaling_factors/calculate_scaling_factors runs once globally
                // (it aggregates every plate via plate_measurements...collect() below),
                // not once per plate like index_unscaled - so finalize_fp.fingerprint
                // (one emission per plate) needs collapsing into a single deterministic
                // value first, reusing the same per-well fingerprints computed above.
                finalize_fp_global = finalize_fp.fingerprint.collect().map{ it.flatten().sort() }
                    
                estimate_scaling_factors(
                    finalize_fp_global,
                    plate_measurements.map{row -> row[1]}.collect(),
                    sc_autoscale,
                    sc_manualscale,
                    sc_scale_slope,
                    sc_scale_bias,
                    sc_control_list,
                    sc_channel_map
                )
                scaling_file = estimate_scaling_factors.out.scaling_file
                slope_file = estimate_scaling_factors.out.slope_file
                bias_file = estimate_scaling_factors.out.bias_file
            } else if (sc_manualscale != null){
                // Manual scaling
                scaling_file = file(sc_manualscale)
                slope_file = (sc_scale_slope != null) ? file(sc_scale_slope) : file("NO_SLOPE")
                bias_file = (sc_scale_bias != null)  ? file(sc_scale_bias)  : file("NO_BIAS")
            } 
            
            // Rescale the existing cached processed images
            rescale_out = rescale(finalize_out.processed_output, scaling_file, slope_file, bias_file)

            // rescale always writes to rr__processed_images/scaled - index it separately
            // from the unscaled output indexed above. See buildPlateFingerprint above
            // for why previous_completed carries a content fingerprint instead of a
            // throwaway barrier value.
            rescale_fp = buildPlateFingerprint(rescale_out.processed_output)
            index_scaled(rescale_fp.fingerprint,
                            "rr__processed_images/scaled",
                            file(rn_publish_dir + "/rr__processed_images/scaled"),
                            rescale_fp.plate)
                            
            images_out = rescale_out.processed_output
        } else {
            // Already wrote the final scaled images directly above - no unscaled artifact
            // exists to measure or estimate factors from, and no rescale is needed
            images_out = finalize_out.processed_output
        } 
        
        //------------------------------------------------------------------------
        // Block 3: paths converge again, and the output is used to make cellcrops
        //------------------------------------------------------------------------
        
        // Combined channel of well, images, masks - used by downstream steps that need both
        finalize_images_and_masks = images_out.join(finalize_out.mask_output, by: 0)

        if (rn_make_cellcrops) {
            if (rn_manifest_registration != null) {
                // If we are registering cycles, we need to update the manifest to reflect the new channel indices
                // This is used to calculate the cycle correlations
                // Construct the channel to update the indices of <plate>:<old_channel> <registered_channel>
                // Used .unique() before, but replaced with custom logic so it doesn't need to wait for the channel to complete
                def seen = []
                channel_map = finalize_out.channel_indices
                .map { item ->
                    if( !seen.contains(item) ) {
                        seen << item
                        return item
                    }
                    return null
                }
                .filter { it != null }
                .flatMap{ manifest_path -> file(manifest_path)
                .splitCsv(header:["ref_plate", "plate", "cycle", "channel", "name", "orig_channel", "orig_name"], sep:"\t")
                }
                .map{ row -> tuple(row.plate + ":" + row.orig_channel, row.channel) }

                // Create a new registration channel where the channel has been updated to the post-registration channel index
                manifest_registration_updated = manifest_registration
                .flatMap {row ->
                    def result = []

                    // Add the ref plate
                    result <<  tuple(row.ref_plate + ":" + row.ref_channel, groupKey(row.ref_plate, row.qry_channels.size()+1), row.ref_channel, null, null)

                    // Add the query plates
                    row.qry_plates.eachWithIndex{
                        val, idx ->
                        result <<  tuple(val + ":" + row.qry_channels[idx], groupKey(row.ref_plate, row.qry_channels.size()+1), row.ref_channel, val, row.qry_channels[idx])
                    }
                    return result
                }
                .combine(channel_map, by: 0)
                .map{ row -> tuple(row[1], row[2], row[3], row[5]) }
                .groupTuple(by: 0)
                .map{ row -> RegistrationRecord.fromVars(row[0].getGroupTarget(), row[1][0], row[2].findAll(), row[3].findAll{ idx -> row[2][row[3].indexOf(idx)] != null}) } // findall removes the null (refplates) so they are not treated as qry

                // Create the cellcrop input channel with the updated registration record that has the registered channel indices to correlate
                cellcrop_in = finalize_images_and_masks
                .map{ row -> tuple(row[0].plate, row[0], row[1], row[2]) }
                .combine(manifest_registration_updated.map{ row -> tuple(row.ref_plate, row) }, by: 0)
                .map{ row -> tuple(row[1], row[4], row[2], row[3]) }
            } else {
                // No registration, just pass through, setting registration to null
                cellcrop_in = finalize_images_and_masks.map{ row -> tuple(row[0], null, row[1], row[2]) }
            }

            // Create cellcrops
            cellcrop_out = cellcrops(cellcrop_in)

            // Index cellcrops - see buildCellcropsFingerprint above for why this isn't
            // cellcrop_out.h5.last() (an arbitrary last-completed well, order-unstable
            // and blind to in-place content changes several directories deep)
            index_cellcrops(buildCellcropsFingerprint(cellcrop_out.h5), file(rn_publish_dir + "/rr__cellcrops"))
        }

    emit:
        finalize_images_and_masks = finalize_images_and_masks
}
