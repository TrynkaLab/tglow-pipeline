#!/usr/bin/env nextflow

include { finalize; rescale; index_cellcrops; index_cellcrops_csv; cellcrops } from '../processes/finalize.nf'
include { measure_intensity; stage_as_plate } from '../processes/intensity.nf'
include { index_images as index_unscaled; index_images as index_scaled } from '../processes/staging.nf'
include { estimate_scaling_factors } from './scaling_factors.nf'

// Finalizes and caches images for feature extraction (rn_cache_images).
//
// Normal flow (sc_scale_in_finalize=false): finalize writes unscaled images to
// processed_images/unscaled, intensities are always measured on those (even if no scaling is
// requested, to support a future QC step), then - if scaling is requested - factors are
// estimated from those same unscaled images and a separate rescale pass produces
// processed_images/scaled.
//
// Fast path (sc_scale_in_finalize=true, requires sc_manualscale - enforced in checks.nf):
// scaling factors are already known upfront, so finalize applies them directly and writes
// straight to processed_images/scaled, skipping the unscaled artifact, measurement, and the
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
        index_unscaled(finalize_out.processed_output.last(),
                        "processed_images/${finalize_subdir}",
                        file(rn_publish_dir + "/processed_images/${finalize_subdir}"),
                        finalize_out.processed_output.map{row -> row[0].plate}.unique())

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

                estimate_scaling_factors(
                    finalize_out.processed_output.last(),
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

            // rescale always writes to processed_images/scaled - index it separately
            // from the unscaled output indexed above
            index_scaled(rescale_out.processed_output.last(),
                            "processed_images/scaled",
                            file(rn_publish_dir + "/processed_images/scaled"),
                            rescale_out.processed_output.map{row -> row[0].plate}.unique())
                            
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
                .map{ row -> new RegistrationRecord(row[0].getGroupTarget(), row[1][0], row[2].findAll(), row[3].findAll{ idx -> row[2][row[3].indexOf(idx)] != null}) } // findall removes the null (refplates) so they are not treated as qry

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

            // Index cellcrops
            index_cellcrops(cellcrop_out.h5.last(), file(rn_publish_dir + "/cellcrops"))
        }

    emit:
        finalize_images_and_masks = finalize_images_and_masks
}
