#!/usr/bin/env nextflow

include { finalize; rescale; index_cellcrops; cellcrops } from '../processes/finalize.nf'
include { measure_intensity } from '../processes/scaling.nf'
include { index_images as index_unscaled; index_images as index_scaled } from '../processes/staging.nf'
include { estimate_scaling_factors } from './scaling_factors.nf'

import RegistrationRecord

// Finalizes and caches images for feature extraction (rn_cache_images).
//
// Normal flow (rn_scale_in_finalize=false): finalize writes unscaled images to
// processed_images/unscaled, intensities are always measured on those (even if no scaling is
// requested, to support a future QC step), then - if scaling is requested - factors are
// estimated from those same unscaled images and a separate rescale pass produces
// processed_images/scaled.
//
// Fast path (rn_scale_in_finalize=true, requires rn_manualscale - enforced in checks.nf):
// scaling factors are already known upfront, so finalize applies them directly and writes
// straight to processed_images/scaled, skipping the unscaled artifact, measurement, and the
// separate rescale pass entirely.
workflow finalize_images {

    take:
        finalize_in
        image_dir_file
        flatfield_out
        manifest_registration
        rn_scale_in_finalize
        rn_autoscale
        rn_manualscale
        rn_scale_slope
        rn_scale_bias
        rn_make_cellcrops
        rn_manifest_registration
        rn_publish_dir
        // Ingredients for estimating scaling factors from the unscaled finalized images
        // (only used when rn_scale_in_finalize is false)
        manifest
        blacklist_file
        control_file
        plates
        manifest_registration_file
        cellpose_out
        rn_control_list

    main:

        // Determine if scaling needs to be run at all
        run_scaling = rn_manualscale != null || rn_autoscale

        //------------------------------------------------------------------------
        // Block 1: finalize unscaled OR direct manual scaled images wihtout caching them in between
        //------------------------------------------------------------------------
        if (rn_scale_in_finalize && run_scaling) {
            finalize_scaling_in = Channel.value(file(rn_manualscale))
            finalize_slope_in   = (rn_scale_slope != null) ? Channel.value(file(rn_scale_slope)) : Channel.value(file("NO_SLOPE"))
            finalize_bias_in    = (rn_scale_bias != null)  ? Channel.value(file(rn_scale_bias))  : Channel.value(file("NO_BIAS"))
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

        finalize_subdir = rn_scale_in_finalize ? "scaled" : "unscaled"

        // Index whatever finalize actually produced (index_unscaled despite the name -
        // when rn_scale_in_finalize writes directly to `scaled`, this call indexes that
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
        if (!rn_scale_in_finalize && run_scaling) {
            if (rn_autoscale) {
                // Auto scaling
                estimate_scaling_factors(
                    manifest,
                    manifest_registration,
                    blacklist_file,
                    control_file,
                    plates,
                    manifest_registration_file,
                    cellpose_out,
                    finalize_out.processed_output.last(),
                    flatfield_out,
                    rn_autoscale,
                    rn_manualscale,
                    rn_manifest_registration,
                    rn_scale_slope,
                    rn_scale_bias,
                    rn_control_list,
                    rn_publish_dir
                )
                scaling_file = estimate_scaling_factors.out.scaling_file
                slope_file = estimate_scaling_factors.out.slope_file
                bias_file = estimate_scaling_factors.out.bias_file
            else if (rn_manualscale != null){
                // Manual scaling
                scaling_file = file(rn_manualscale)
                slope_file = (rn_scale_slope != null) ? file(rn_scale_slope) : file("NO_SLOPE")
                bias_file = (rn_scale_bias != null)  ? file(rn_scale_bias)  : file("NO_BIAS")
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
                .map(row -> tuple(row.plate + ":" + row.orig_channel, row.channel))

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
                .map(row -> tuple(row[1], row[2], row[3], row[5]))
                .groupTuple(by: 0)
                .map(row -> new RegistrationRecord(row[0].getGroupTarget(), row[1][0], row[2].findAll(), row[3].findAll{ idx -> row[2][row[3].indexOf(idx)] != null})) // findall removes the null (refplates) so they are not treated as qry

                // Create the cellcrop input channel with the updated registration record that has the registered channel indices to correlate
                cellcrop_in = finalize_images_and_masks
                .map(row -> tuple(row[0].plate, row[0], row[1], row[2]))
                .combine(manifest_registration_updated.map(row -> tuple(row.ref_plate, row)), by: 0)
                .map(row -> tuple(row[1], row[4], row[2], row[3]))
            } else {
                // No registration, just pass through, setting registration to null
                cellcrop_in = finalize_images_and_masks.map(row -> tuple(row[0], null, row[1], row[2]))
            }

            // Create cellcrops
            cellcrop_out = cellcrops(cellcrop_in)

            // Index cellcrops
            index_cellcrops(cellcrop_out.h5.last(), file(rn_publish_dir + "/cellcrops"))
        }

    emit:
        finalize_images_and_masks = finalize_images_and_masks
}
