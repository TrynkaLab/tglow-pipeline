#!/usr/bin/env nextflow

// Processes
include { convertChannelType; readBlacklist} from '../lib/utils.nf'
include { deconvolute } from '../processes/decon.nf'
include { register } from '../processes/registration.nf'
include { cellpose } from '../processes/segmentation.nf'
include { cellprofiler;  finalize_and_cellprofiler } from '../processes/cellprofiler.nf'
include { finalize; index_cellcrops; cellcrops } from '../processes/finalize.nf'
include { index_imagedir } from '../processes/staging.nf'

// Subworkflows
include { setup } from "../subworkflows/setup.nf"
include { flatfield_estimation } from "../subworkflows/flatfield_estimation.nf"
include { estimate_scaling_factors } from "../subworkflows/scaling_factors.nf"

import Well
import ManifestRecord
import RegistrationRecord

// Main workflow
workflow run_pipeline {

    main:
    
        //------------------------------------------------------------
        // Check the input parameters
        //------------------------------------------------------------
        // Sanity check input
        if (params.cpr_run) {
            if (params.rn_max_project || params.rn_hybrid) {
                if (params.cpr_pipeline_2d == null) {
                    error("Running in hybrid or 2d mode with cpr_run=true, but cpr_pipeline_2d is not set")
                }
                
                if (!file(params.cpr_pipeline_2d).exists()) {
                    error("Specified cellprofiler pipeline is not accessible")
                }
                
            } else {
                if (params.cpr_pipeline_3d == null) {
                    error("Running in 3d mode with cpr_run=true, but cpr_pipeline_3d is not set")
                }
                
                if (!file(params.cpr_pipeline_3d).exists()) {
                    error("Specified cellprofiler pipeline is not accessible")
                }
                
            }
        }
        
        // Check cellcrops
        if (params.rn_make_cellcrops && !params.rn_cache_images) {
            error("rn_cache_images must be true when rn_make_cellcrops is true")
        }

        //------------------------------------------------------------
        // Run setup
        //------------------------------------------------------------
        // Run the setup, parsing manifests and setting up channels
        setup(params)
        
        // Channels with file contents
        manifest = setup.out.manifest
        plates = setup.out.plates
        manifest_registration = setup.out.manifest_registration
        
        // Channels with file objects
        manifest_registration_file = setup.out.manifest_registration_file
        blacklist_file = setup.out.blacklist_file
        control_file = setup.out.control_file

        //------------------------------------------------------------
        // Prepare per well input channels
        //------------------------------------------------------------
        // Loop over previously generated manifests assuming stage has been run
        if (params.rn_manifest_well == null) {            
            manifests_in = manifest.map{row -> "${params.rn_image_dir}/" + row.plate + "/manifest.tsv"}
        } else {
            manifests_in = Channel.from(params.rn_manifest_well.split(','))
        }
        
        // Construct the channel on the well level
        well_channel = manifests_in
            .flatMap{ manifest_path -> file(manifest_path)
            .splitCsv(header:["well", "row", "col", "plate"], sep:"\t")}
            .map( row -> new Well(well: row.well, row: row.row, col: row.col, plate: row.plate) )

                
        // Filter blacklist. Blacklist read into arrat of <plate>:<well>
        if (params.rn_blacklist != null) {
            blacklist = readBlacklist(params.rn_blacklist)
            well_channel = well_channel.filter(row -> {row.key !in blacklist})
        }
            
        // Filter to specific wells, usefull for testing
        if (params.rn_wells != null) {
            wells = params.rn_wells.split(",")
            well_channel = well_channel.filter(row -> {row.well in wells})
        }                  

        // Add the plate for easier combining later
        well_channel = well_channel.map{ row -> tuple(row.plate, row) }
        
        //------------------------------------------------------------
        // Run flatfield estimation
        //------------------------------------------------------------
        flatfield_estimation(
            manifest,
            manifest_registration,
            blacklist_file,
            params.rn_manifest,
            params.rn_manifest_registration,
            params.bp_run,
            params.bp_channels,
            params.bp_global_flatfield,
            plates
        )
        
        flatfield_out = flatfield_estimation.out.flatfield_out
        
        //------------------------------------------------------------
        // Deconvolute
        //------------------------------------------------------------
        def image_dir_file = file(params.rn_image_dir)
        if (params.dc_run) {
            
            decon_in = well_channel
            .combine(manifest.map{row -> tuple(row.plate, row)}, by: 0)
            .filter(row -> {row[2].dc_channels != "none"})
            .map{ row -> tuple(
                row[1], // Well
                row[2], // Manifest
                row[2].dc_channels.collect(it -> {it.toString() + "=" + new File(row[2].dc_psfs[it]).getName()}).sort().join(" "), // String for channel file pairs
                row[2].dc_psfs.collect(it -> file(it)), // The psf files
                file(params.rn_image_dir + "/" + row[1].relpath) // The image files
            )}
            
            // A channel of Well, ManifestRecord, path to images
            decon_out = deconvolute(decon_in)

            cellpose_in = decon_out.filter(row -> {row[1].cp_cell_channel != "none"})

            // The image dir is now the decon output
            image_dir_file = file(params.rn_decon_dir)
        } else {
            // Channel of Well, ManifestRecord, path to images
            cellpose_in = well_channel
            .combine(manifest.map{row -> tuple(row.plate, row)}, by: 0)
            .filter(row -> {row[2].cp_cell_channel != "none"})
            .map{row -> tuple(row[1], row[2], file(params.rn_image_dir + "/" + row[1].relpath))} 
        }
        
        //------------------------------------------------------------
        // Registeration
        //------------------------------------------------------------
        if (params.rn_manifest_registration != null) {
            // Well, RegistrationRecord, path to images
            registration_in = well_channel
            .combine(manifest_registration.map{row -> tuple(row.ref_plate, row)}, by:0)
            .map{ row -> tuple(row[1], row[2], file(params.rn_image_dir))} 

            // Run registration
            registration_out = register(registration_in)
                          
            // Filter cellpose channel to run reference plates only 
            // Channel of Well, ManifestRecord, path to images
            cellpose_in = cellpose_in.map(row -> tuple(row[0].plate, row))
            .combine(manifest_registration.map{row -> tuple(row.ref_plate, row)}, by: 0)
            .map{ row -> tuple(row[1])}
             
        } else {
            registration_out = null
        }
        
        
        //------------------------------------------------------------
        // Cellpose
        //------------------------------------------------------------
        // Run cellpose
        if (params.cp_run) {
            cellpose_out = cellpose(cellpose_in)
        } else {
            cellpose_out = Channel.empty()
        }
        
        //------------------------------------------------------------
        // Estimate scaling factors
        //------------------------------------------------------------
        estimate_scaling_factors(
            manifest,
            manifest_registration,
            blacklist_file,
            control_file,
            plates,
            manifest_registration_file,
            cellpose_out,
            flatfield_out,
            params.rn_autoscale,
            params.rn_manualscale,
            params.rn_manifest_registration,
            params.rn_scale_slope,
            params.rn_scale_bias
        )
        
        scaling_file = estimate_scaling_factors.out.scaling_file
        slope_file = estimate_scaling_factors.out.slope_file
        bias_file = estimate_scaling_factors.out.bias_file

        //------------------------------------------------------------
        // Finalize
        //------------------------------------------------------------

        // This ensures that all the decons are done before finalize starts
        if (params.rn_manifest_registration != null && params.dc_run) {

            // qry plate, ref plate (used as grouping key)
            // Setting the size of the grouping key allows it to emit items as soon as the decons are done
            qry_ref =  manifest_registration.flatMap{
                pair ->
                def key = pair.ref_plate
                def values= pair.qry_plates
                def result = [[key, groupKey(key, values.size()+1)]]
                values.each { val ->
                    result << [val, groupKey(key, values.size()+1)]
                }
                return result
            }

            // Add the reference plate to the decon output
            // Generates a channel, ref_plate:well, cyle plates [decon_c1, decon_c2, ...]
            // This is not really used other that for ensuring that both decons are done before
            // I tried to specify the files instead of the plates, but this is hard to stage
            // in a plate/row/col structure for finalize as the file directories have the same name
            image_input = decon_out
                .map{row -> tuple(row[0].plate, row[0], row[1], row[2])}
                .combine(qry_ref, by: 0)
                .map(row -> tuple(groupKey(row[4] + ":" + row[1].well, row[4].getGroupSize()), row[0], row[1], row[2], row[3])) // ref plate, plate, well, manifest, path
                .groupTuple(by: 0)
                .map(row -> tuple(row[0].getGroupTarget(), row[1]))
        } else {
            // Generates a channel, ref_plate:well, [plates]
            image_input = cellpose_in.map(row -> tuple(row[0].key, row[1].plate)) // key, cycle plates
        }
        
        //--------------------------------------------------------------------
        // Add the cell and nucleus masks for each well
        // key, Well, ManifestRecord, cell masks, nucl masks,
        finalize_in = cellpose.out
            .map(row -> tuple(row[0].key, row[0], row[1], row[2], row[3]))
            .combine(image_input, by:0)
        
        //--------------------------------------------------------------------
        // Add registration
        if (registration_out != null) {
            // re-key output
            // key, merge plates, path
            registration_out = registration_out.map{row -> tuple(row[0].key, row[1].qry_plates, row[2])} 
            finalize_in = finalize_in.join(registration_out, by: 0)                
        } else {
            finalize_in = cellpose_out.map{row -> tuple(
                row[0].key, // key,
                row[0], // Well
                row[1], // ManifestRecord
                row[2], // cell masks,
                row[3], // nucl masks,
                null, // decon plates TODO: figure out what this does, it doesn't seem to be used in finalize or cellprofiler
                null, // merge plates
                file('NO_REGISTRATION')   // registration path
            )}
        }
        
        //--------------------------------------------------------------------
        // Cache the final images for feature extraction
        if (params.rn_cache_images) {
                        
            finalize_out = finalize(finalize_in,
                                    image_dir_file,
                                    flatfield_out,
                                    scaling_file,
                                    slope_file,
                                    bias_file)
                                    
            // Create the plate manfiests once finalize is done
            index_imagedir(finalize_out.processed_output.last(),
                            "processed_images",
                            file(params.rn_publish_dir + "/processed_images"),
                            finalize_out.processed_output.map{row -> row[0].plate}.unique())
                        
            if (params.rn_make_cellcrops) {
                 
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

                cellcrop_in = finalize_out.processed_output
                .map(row -> tuple(row[0].plate, row[0], row[1]))
                .combine(manifest_registration_updated.map(row -> tuple(row.ref_plate, row)), by: 0)
                .map(row -> tuple(row[1], row[3], row[2]))
                
                // Create cellcrops
                cellcrop_out = cellcrops(cellcrop_in)
                
                // Index cellcrops
                cellcrop_out.h5.last().ifPresent { h5 ->
                    index_cellcrops(h5, file(params.rn_publish_dir + "/cellcrops"))
                }  
            }            
        }
        
        //------------------------------------------------------------
        //                       Cellprofiler
        //------------------------------------------------------------
        // Run cellprofiler        
        if (params.cpr_run && params.rn_cache_images) {
            // Run cellprofiler on cached images
            cellprofiler_out = cellprofiler(finalize_out.processed_output)
        } else if (params.cpr_run) {
            // This does not cache images
            cellprofiler_out = finalize_and_cellprofiler(finalize_in,
                                                        image_dir_file,
                                                        flatfield_out,
                                                        scaling_file,
                                                        slope_file, 
                                                        bias_file)
        } 
        
         
}
