#!/usr/bin/env nextflow

include { convertChannelType; parseManifestFlatfields; readBlacklist} from '../lib/utils.nf'
include { deconvolute } from '../processes/decon.nf'
include { estimate_flatfield; stage_global_flatfield } from '../processes/flatfields.nf'
include { register } from '../processes/registration.nf'
include { calculate_scaling_factors; calculate_plate_offsets } from '../processes/scaling.nf'
include { cellpose } from '../processes/segmentation.nf'
include { cellprofiler; finalize; finalize_and_cellprofiler } from '../processes/features.nf'
include { index_imagedir } from '../processes/staging.nf'


// Main workflow
workflow run_pipeline {

    main:
        //------------------------------------------------------------
        // Defaults & sanity checks
        //------------------------------------------------------------
        if (params.rn_publish_dir == null) {
            error "rn_publish_dir file parameter is required: --rn_publish_dir"
        }
        
        if (params.rn_manifest == null) {
            error "rn_manifest file parameter is required: --rn_manifest"
        }
        
        if (params.rn_cache_images & !(params.rn_max_project | params.rn_hybrid)) {
            log.warn "Caching images in 3d mode can take up a lot of space, are you sure this is what you want?"
        }
        
        if (!params.rn_autoscale & params.rn_control_list != null) {
            error "Provided --rn_control_list but --rn_autoscale false. Either drop --rn_control_list or set --rn_autoscale true"
        }
        
        if ((params.rn_manualscale != null) & params.rn_autoscale) {
            log.warn "Both rn_autoscale and rn_manualscale are provided, rn_manualscale will be overridden"
        }


        // Set runtime defaults, these are overidden when specified on commandline
        params.rn_image_dir = params.rn_publish_dir + "/images"
        params.rn_decon_dir = params.rn_publish_dir + "/decon"
                
        //------------------------------------------------------------
        // Read manifests & input files
        //------------------------------------------------------------
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(
            row.plate,
            row.index_xml,
            (row.channels == null) ? "none" : tuple(row.channels.split(',')), 
            (row.bp_channels == null) ? "none" : tuple(row.bp_channels.split(',')),
            row.cp_nucl_channel,
            row.cp_cell_channel,
            row.dc_channels,
            row.dc_psfs,
            (row.mask_channels == null) ? "none" : tuple(row.mask_channels.split(',')))}
        
        // Build a value channel for the plates in the manifest
        // '<plate_1> <plate_2> <plate_N>'
        plates_channel = manifest.map{row -> row[0]}.collect().map { it.join(' ') }
        
        //------------------------------------------------------------------------
        // Registration manifest
        if (params.rn_manifest_registration != null) {
            manifest_registration = Channel
            .fromPath(params.rn_manifest_registration)
            .splitCsv(header:true, sep:"\t")
            .map{row -> tuple(
                row.reference_plate,
                row.reference_channel,
                row.query_plates,
                row.query_channels
            )}               
        } 
    
        //------------------------------------------------------------------------
        // Blacklist channel, if missing just an empty channel
        if (params.rn_blacklist == null) {
            log.info("No blacklist provided")
            blacklist_channel = Channel.value(file('NO_BLACKLIST'))
        } else {
            blacklist_channel = Channel.value(file(params.rn_blacklist))
        }
        
        //------------------------------------------------------------------------
        // Control list channel, if missing just an empty channel
        if (params.rn_control_list == null) {
            log.info("No controlist provided")
            control_list_channel = Channel.value(file('NO_CONTROL_LIST'))
        } else {
            control_list_channel = Channel.value(file(params.rn_control_list))
        }
    
        //------------------------------------------------------------------------
        // Registration manifest, if missing just an empty channel
        if (params.rn_manifest_registration == null) {
            log.info("No registration provided")
            manifest_registration_channel = Channel.value(file('NO_REGISTRATION_MANIFEST'))
        } else {
            manifest_registration_channel = Channel.value(file(params.rn_manifest_registration))
        }
    
        //------------------------------------------------------------
        // Run estimate_flatfield
        //------------------------------------------------------------
        // If there is no global overide on flatfield channels, get them from manifest
        def run_flatfield = false
        if (params.bp_channels == null) {
            // Manually parse the manifest to also grab the first record which is hard to do with channels
            flatfield_channels = parseManifestFlatfields(params.rn_manifest)
            
            if (flatfield_channels[0].size() > 0) {run_flatfield = true}
            
            // These have the channel already zero indexed
            flatfield_in = Channel.from(flatfield_channels[0])
            flatfield_in_global = Channel.from(flatfield_channels[1])
                            
            if (params.rn_manifest_registration) {
                log.info("Estimating one global flatfield per cycle")
                
                // Create a flat [<plate> <cycle>] channel
                plate_cycle = manifest_registration
                .map{row -> [row[0], *row[2].split(",")]}
                .flatMap { row -> row.withIndex().collect { item, index ->[item, index] }}
                
                // Create a channel per cycle <ref_plate> <cycle> [<plate1> <plateN>]
                per_cycle = plate_cycle
                .groupTuple(by:1)
                .map{row -> tuple(row[0][0], row[1], row[0])}
                
                // Subtract one from the channel here, as we combine the manifest which is one indexed
                flatfield_in_global = per_cycle
                .combine(manifest, by:0)
                .flatMap{row -> row[5].collect{ channel -> tuple(row[1] + ":" + (channel.toInteger()-1).toString(), row[1], row[2][0], row[2],  channel.toInteger() - 1, row[3])}} // key, cycle, plate, plate(s), channel, index xml
                
                // These have already been converted to zero indexed
                flatfield_in = flatfield_in
                .combine(plate_cycle, by:0)
                .map{row -> tuple(row[3] + ":" + row[1], row[3], row[0], [row[0]],row[1], row[2])} // key, cycle, plate, plate(s), channel, index xml
                                
            } else {
                
                flatfield_in_global = flatfield_in_global
                .combine(plates_channel)
                .map{row -> tuple( "0:" + row[1], 0, row[0], row[3].split(" "), row[1], row[2])} // key, cycle, plates, channel, index_xml
                
                flatfield_in = flatfield_in
                .combine(plates_channel)
                .map{row -> tuple("0:" + row[1], 0, row[0], row[3].split(" "), row[1], row[2])} // key, cycle, plates, channel, index_xml
            }
        } else {
            flatfield_in = Channel.from(params.bp_channels)
        }
        
        if (run_flatfield) {
            log.info("Running flatfield estimation")
            if (params.bp_global_flatfield) {

                // Runs only one model per plate
                global_flatfield = estimate_flatfield(flatfield_in_global, blacklist_channel).flatfield_out
                
                // Now need to make a channel like flatfield_in, but replacing the path with global_flatfield
                // Combine on the channel
                global_flatfield_in = flatfield_in
                .combine(global_flatfield, by: 0)
                .map{ row -> tuple(
                    row[0], // key
                    row[1], // cycle
                    row[2], // plate
                    row[3], // plates
                    row[4], // channel
                    (row[10] instanceof List) ? row[10].findAll{it.fileName.name.startsWith("global_refplate")} : row[10] // reference path
                )}        
                
                // I THINK THE BELOW IS FIXED NOW
                // This was the old way before, I don't remember exactly why I did this, but its giving issues when not registering
                // I think this was to avoid cases where there has been a previous pipeline run
                // Need to make this check conditional on the case there is a multi return
                //  row[10].findAll{it.fileName.name.startsWith("global_refplate")} // reference path
                //global_flatfield_in.view()        
                
                flatfield_out = stage_global_flatfield(global_flatfield_in).flatfield_out
            } else {
                log.info("Estimating one flatfield per channel")
                flatfield_out = estimate_flatfield(flatfield_in, blacklist_channel).flatfield_out
            }
            
            // Concat to plate string format // plate : channel = path
            flatfield_out_string = flatfield_out.map{ row -> (row[2] + "_ch" + row[4] + "=" + row[5])}.collect().map{it.join(" ")}
        } else {
            flatfield_out_string = Channel.value(false)
            log.info("No manifest entries or --bp_channels provided, so skipping flatfield estimation")
        }
         
        //------------------------------------------------------------
        // Prepare per well input channels
        //------------------------------------------------------------
        // Loop over previously generated manifests assuming stage has been run
        if (params.rn_manifest_well == null) {            
            manifests_in = manifest.map{row -> "${params.rn_image_dir}/" + row[0] + "/manifest.tsv"}
        } else {
            manifests_in = Channel.from(params.rn_manifest_well.split(','))
        }
        
        // Construct the channel on the well level
        well_channel = manifests_in
            .flatMap{ manifest_path -> file(manifest_path)
            .splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }
                                 
        // Filter blacklist. Blacklist read into arrat of <plate>:<well>
        if (params.rn_blacklist != null) {
            blacklist = readBlacklist(params.rn_blacklist)
            well_channel = well_channel.filter(row -> {(row.plate + ":" + row.well) !in blacklist})
        }
        
        // Filter to specific wells, usefull for testing
        if (params.rn_wells != null) {
            wells = params.rn_wells.split(",")
            well_channel = well_channel.filter(row -> {row.well in wells})
            log.info("Selecting " + wells.size() + " wells: " + wells)
        }                  
    
        // Re-order the well channel for later merging
        well_in = well_channel.map{ row -> tuple(row.plate, row.well, row.row, row.col)}
                
        //------------------------------------------------------------
        // Deconvolute
        //------------------------------------------------------------
        if (params.dc_run) {
            decon_in = well_in.combine(manifest, by: 0)
            .filter(row -> {row[9] != "none"})
            .map{ row -> tuple(
                row[0], row[1], row[2], row[3], // plate, well, row, col
                convertChannelType(row[7]), // nucl_channel
                convertChannelType(row[8]), // cell_channel
                row[9].split(',').collect{it -> (it.toInteger() -1).toString() + "=" + new File(row[10].split(',')[it.toInteger() -1]).getName()}.join(" "),  // dc_channels
                row[10].split(',').collect(it -> (file(it))) //dc_psfs
            )}      
            
            // Deconvolute
            cellpose_in = deconvolute(decon_in)
            .filter(row -> {row[5] != "none"})
            .map{ row -> tuple(
                row[0], row[1], row[2], row[3], // plate, well, row, col
                row[4], // nucl_channel
                row[5] // cell_channel
            )}            
        } else {
            // Append the things from the manifest needed for cellpose
            cellpose_in = well_in.combine(manifest, by: 0)
            .filter(row -> {row[8] != "none"})
            .map{ row -> tuple(
                row[0], row[1], row[2], row[3], // plate, well, row, col
                convertChannelType(row[7]), // nucl_channel
                convertChannelType(row[8]) // cell_channel
            )}
        }
           
        //------------------------------------------------------------
        // Registeration
        //------------------------------------------------------------
        if (params.rn_manifest_registration != null) {
            // Append the plate info from the manifest if it exists
            // Use a join so the plates to register are discarded
            // Use -1 as the channels are 0 indexed
            registration_in = well_in
            .combine(manifest_registration, by:0)
            .map{ row -> tuple(
                row[0], row[1], row[2], row[3], // plate, well, row, col
                row[4].toInteger()-1, // reference_channel
                row[5].split(',').join(" "), // query_plates
                row[6].split(',').collect{it -> (it.toInteger() -1).toString()}.join(" ")  // query_channels
            )} 

            // Run registration
            registration_out = register(registration_in)
                          
            // Filter cellpose channel to run reference plates only 
            cellpose_in = cellpose_in
            .combine(manifest_registration, by: 0)
            .map{ row -> tuple(
                row[0], row[1], row[2], row[3], // plate, well, row, col
                row[4], // nucl_channel
                row[5]  // cell_channel
            )}                
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
        // Scaling
        //------------------------------------------------------------
        // Get scaling string
        if (params.rn_manualscale != null) {
            scaling_channel = Channel
            .fromPath(params.rn_manualscale)
            .map { file ->  file.text}
            .first() // Makes sure it is a value channel
        } else {
            scaling_channel = Channel.value("none")
        }
        
        // When all deconvelution is done, or all data is staged, calculate the scaling factors
        // which map the images onto the full dynamic range of the 16 bit uint
        // Alternatively, if control scaling is provided, wait for cellpose and run the controls.
        if (params.rn_autoscale) {
            // Channel <plate> <mask_channel1 mask_channel2 mask_channelN>
            scaling_in = manifest.map{
                row -> tuple(row[0],
                (row[8][0] == "none") ? "none" : row[8].collect{it -> (row[0] + "=" + (it.toInteger() -1).toString())}.join(" "),
                null)
            }
            
            // Filter to reference plates only as we want scaling factors to be in final dimension
            if (params.rn_manifest_registration != null){
                scaling_in = scaling_in.combine(manifest_registration, by: 0).map{row -> tuple(row[0], row[1], row[4].split(',').join(" "))}
            }
            
            if (params.rn_control_list) {
                // Start only if registration & cellpose is done
                control_dir = calculate_plate_offsets(
                    cellpose_out.last(), // ensures this only starts when cellpose is done
                    scaling_in,
                    Channel.value(file("${params.rn_publish_dir}/registration")),
                    Channel.value(file("${params.rn_publish_dir}/masks")),
                    flatfield_out_string,
                    blacklist_channel,
                    control_list_channel
                ).plate_offset
            } else {
                control_dir = Channel.value(file('NO_CONTROL_DIR'))
            }
            
            // Final scaling
            demultiplex_channelstring = scaling_in.map{row -> row[1]}.collect().map({it.unique().join(" ")})

            // Start running only if all the plate offsets have been calculated
            scaling_channel = calculate_scaling_factors(cellpose_in.last(),
                blacklist_channel,
                plates_channel,
                manifest_registration_channel,
                control_dir.collect(), // ensures this only runs when control dir is done
                demultiplex_channelstring).scaling_factors.map { file -> 
                    file.text
            }
        }


        //------------------------------------------------------------
        // Finalize
        //------------------------------------------------------------
        // TODO: BUG! Need to make sure the merge plate decon is done as well BEFORE adding into channel
                        
        // Start with cellpose output
        // re-key channels
        cellpose_out = cellpose_out.map{row -> tuple(
                row[0] + ":" + row[1], // key
                row[0], // plate
                row[1], // well
                row[2], // row
                row[3], // col
                row[4], // cell masks
                row[5]  // nucl masks
        )}

        //--------------------------------------------------------------------
        // Add registration
        if (registration_out != null) {
            // re-key output
            registration_out = registration_out.map{row -> tuple(row[0] + ":" + row[1], row[4], row[5])} // key, merge plates, path
            // merge
            finalize_in = cellpose_out.join(registration_out, by: 0)                
        } else {
            finalize_in = cellpose_out.map{row -> tuple(
                row[0], // key
                row[1], // plate
                row[2], // well
                row[3], // row
                row[4], // col
                row[5], // cell masks,
                row[6], // nucl masks,
                null,   // merge plates
                file('NO_REGISTRATION')   // registration path
            )}
        }
                    
        //--------------------------------------------------------------------
        // for in hybrid mode
        if (params.rn_hybrid) {
            remapped_manifest = manifest.map{ row -> tuple(
                row[8], // scale channels
                row[0] // plate
            )}

            finalize_in = finalize_in.combine(remapped_manifest, by: 1).map{row -> tuple(
                row[0], // plate
                row[1], // key
                row[2], // well
                row[3], // row
                row[4], // col
                row[5], // cell masks
                row[6], // nucl masks
                row[7], // merge plates
                row[8], // registration path
                (row[9] == "none") ? "none" : row[9].collect{it -> row[0] + "=" + (it.toInteger() -1).toString()}.join(" ") // mask channels   
            )}
            
        } else {
            finalize_in = finalize_in.map{row -> tuple(
                row[1], // plate
                row[0], // key
                row[2], // well
                row[3], // row
                row[4], // col
                row[5], // cell masks
                row[6], // nucl masks
                row[7], // merge plates
                row[8], // registration path
                null // mask channels
            )}
        }
        
        // Cache the final images for feature extraction
        if (params.rn_cache_images) {
            finalize_out = finalize(finalize_in, flatfield_out_string, scaling_channel)[0]
            
            // Create the plate manfiests once finalize is done
            index_imagedir("processed_images", finalize_out.last(), finalize_out.map{row -> row[0]}.unique())
        }
        
        //------------------------------------------------------------
        //                       Cellprofiler
        //------------------------------------------------------------
        // Run cellprofiler
        if (params.cpr_run & params.rn_cache_images) {
            // Run cellprofiler on cached images
            cellprofiler_out = cellprofiler(finalize_out)
        } else if (params.cpr_run) {
            // This does not cache images
            cellprofiler_out = finalize_and_cellprofiler(finalize_in, flatfield_out_string, scaling_channel)
        }      
}




        // DONE: nextflowify this by remaping this from manifest channel
        // This is the old way, cleanup later
        //if (params.rn_manualscale & !params.rn_autoscale) {
            // def csvFile = new File(params.rn_manifest)
            // def csvData = csvFile.readLines()
            // def scaling_string = null

            // def plate_scale = []
            // for (int i = 1; i < csvData.size(); i++) {
            //     def curLine = csvData[i].split('\t')
                
            //     if ((curLine[9] != "none") & (curLine[9] != null)) {
            //         def curChannels = curLine[2].split(',')
            //         for (channel in curChannels) {
            //             curChannel = channel.toInteger()-1
            //             plate_scale << curLine[0] + "_ch" + curChannel + "=" + curLine[9].split(",")[curChannel]
            //         }
            //     }
            // }   


            // TODO: nextflowify this by remaping this from manifest channel
            // code that reads paths available manfiests from previous stage
            //manifests_in = Channel.fromPath("${params.rn_image_dir}/*/manifest.tsv")
            
            // Read which plates should be run, in the case the manifest has been edited
            // after running stage to run not all plates
            // plates = []
            // new File(params.rn_manifest).eachLine{string, index -> 
            //     if (index > 1) {
            //         plates.add(string.split("\t")[0])
            //     }
            // }
            
            // manifest_paths = []
            
            // for (plate in plates) {
            //     manifest_paths += "${params.rn_image_dir}/" + plate + "/manifest.tsv"
            // }
            
            // log.info("Considering plate level manifests: " + manifest_paths)
            // manifests_in = Channel.from(manifest_paths)
            
            
            
            
            // Rather then having it in the tuple, as the flatfields are constant for all wells
            // now its supplied as a value channel, which saves a lot of code
            //   if (run_flatfield) {
            
            //     // Concat to plate string format
            //     flatfield_tmp = flatfield_out.map{ row -> tuple(
            //         row[0], // plate
            //         row[0] + "_ch" + row[1] + "=" + row[2] // plate : channel = path
            //     )}
    
            //     // Merge channel of same plates into single string 
            //     flatfield_tmp = flatfield_tmp
            //     .groupTuple(by:0)
            //     .map{row -> tuple(
            //         row[1].join(" "), //  plate_ch1=path1 plate_ch2=path2 plate_chX=pathX
            //         row[0] //plate
            //     )}
                
            //     // If there is a registration manifest, the basicpy output needs to be combined
            //     // into a single channel with the reference plate as the key
            //     if (params.rn_manifest_registration) {
                
            //         // Read the manifest and turn it into a (qry, ref) channel
            //         qry_ref = parseRegManfiestAsQryRef(params.rn_manifest_registration)
            //         qry_ref_channel = Channel.from(qry_ref).map{row -> tuple(row[1], row[0])} // ref, qry
                
            //         // Add the reference plate to the basicpy output channel: (plate, bp_string, ref_plate)
            //         flatfield_tmp = flatfield_tmp.combine(qry_ref_channel, by:1)
                    
            //         // Group together by reference plate and concat the bp strings
            //         flatfield_tmp = flatfield_tmp
            //         .groupTuple(by:2)
            //         .map{row -> tuple(
            //             row[1].join(" "), //  plate_ch1=path1 plate_ch2=path2 plate_chX=pathX
            //             row[2] //plate
            //         )}
            //     }

            //     // Append the basicpy models for a plate into that channel              
            //     cellprofiler_in = cellprofiler_in.combine(flatfield_tmp, by: 1)
            // } else {
            //     cellprofiler_in = cellprofiler_in.map{row -> tuple(
            //             row[1], // plate
            //             row[0], // key
            //             row[2], // well
            //             row[3], // row
            //             row[4], // col
            //             row[5], // cell masks
            //             row[6], // nucl masks
            //             row[7], // merge plates
            //             row[8], // registration path
            //             null    // basicpy models       
            //     )}
                
            // }