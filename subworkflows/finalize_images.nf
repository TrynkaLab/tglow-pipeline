#!/usr/bin/env nextflow

include { finalize; rescale; index_cellcrops; cellcrops } from '../processes/finalize.nf'
include { index_images as index_unscaled; index_images as index_scaled } from '../processes/staging.nf'

import RegistrationRecord

// Finalizes and caches images for feature extraction (rn_cache_images), applying scaling
// either directly in `finalize` or via a separate `rescale` pass, and optionally builds
// cellcrops.
workflow finalize_images {

    take:
        finalize_in
        image_dir_file
        flatfield_out
        scaling_file
        slope_file
        bias_file
        manifest_registration
        rn_scale_in_finalize
        rn_autoscale
        rn_manualscale
        rn_make_cellcrops
        rn_manifest_registration
        rn_publish_dir

    main:
    
        // Whether finalize actually ends up writing to processed_images/scaled directly,
        // uniform for the whole run - see subworkflows/finalize.nf for how this is used.
        finalize_writes_scaled = rn_scale_in_finalize && (rn_autoscale || rn_manualscale != null)

        // finalize_writes_scaled mirrors the img_subdir logic inside `finalize` itself -
        // which folder it actually wrote to for this run (passed in by the caller since it's
        // uniform across all wells, only depending on run-level config, not per-well data).
        finalize_subdir = finalize_writes_scaled ? "scaled" : "unscaled"

        // rn_scale_in_finalize: apply scaling directly in finalize (writes straight to
        // processed_images/scaled, skips the separate rescale pass) - trades away the
        // unscaled artifact and finalize/scaling-factor parallelism for less IO on very
        // IO-heavy runs. Default is off: finalize writes to processed_images/unscaled and
        // rescale (if scaling is enabled) produces processed_images/scaled from that.
        finalize_scaling_in = rn_scale_in_finalize ? scaling_file : Channel.value(file("NO_SCALE"))
        finalize_slope_in   = rn_scale_in_finalize ? slope_file   : Channel.value(file("NO_SLOPE"))
        finalize_bias_in    = rn_scale_in_finalize ? bias_file    : Channel.value(file("NO_BIAS"))

        finalize_out = finalize(finalize_in,
                                image_dir_file,
                                flatfield_out,
                                finalize_scaling_in,
                                finalize_slope_in,
                                finalize_bias_in)

        // Index whatever finalize actually produced (index_unscaled despite the name -
        // when rn_scale_in_finalize writes directly to `scaled`, this call indexes that
        // instead; it's just the alias used for "whatever finalize wrote" since a process
        // can only be invoked once per workflow scope)
        index_unscaled(finalize_out.processed_output.last(),
                        "processed_images/${finalize_subdir}",
                        file(rn_publish_dir + "/processed_images/${finalize_subdir}"),
                        finalize_out.processed_output.map{row -> row[0].plate}.unique())

        // Scaling is enabled (scaling_file is real, not the NO_SCALE placeholder) but finalize
        // didn't apply it directly - defer to a separate rescale pass
        if (!rn_scale_in_finalize && scaling_file.name != "NO_SCALE") {
            rescale_out = rescale(finalize_out.processed_output, scaling_file, slope_file, bias_file)
            images_out = rescale_out.processed_output

            // rescale always writes to processed_images/scaled - index it separately
            // from the unscaled output indexed above
            index_scaled(rescale_out.processed_output.last(),
                            "processed_images/scaled",
                            file(rn_publish_dir + "/processed_images/scaled"),
                            rescale_out.processed_output.map{row -> row[0].plate}.unique())
        } else {
            images_out = finalize_out.processed_output
        }

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
