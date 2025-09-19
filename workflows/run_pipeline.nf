#!/usr/bin/env nextflow

// Processes
include { convertChannelType; readBlacklist} from '../lib/utils.nf'
include { deconvolute } from '../processes/decon.nf'
include { register } from '../processes/registration.nf'
include { calculate_scaling_factors; calculate_plate_offsets } from '../processes/scaling.nf'
include { cellpose } from '../processes/segmentation.nf'
include { cellprofiler; finalize; finalize_and_cellprofiler } from '../processes/features.nf'
include { index_imagedir } from '../processes/staging.nf'

// Subworkflows
include { setup } from "../subworkflows/setup.nf"
include { flatfield_estimation } from "../subworkflows/flatfield_estimation.nf"

// Main workflow
workflow run_pipeline {

    main:

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
            params.bp_global_flatfield
        )
        
        flatfield_out = flatfield_estimation.out.flatfield_out
        flatfield_out_string = flatfield_estimation.out.flatfield_out_string

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
            scaling_channel = channel.value(file(params.rn_manualscale))
        } else {
            scaling_channel = Channel.value(file("NO_SCALE"))
        }
        
        // Optional slope for the sigmoid curve
        if (params.rn_scale_slope != null) {    
            slope_channel = channel.value(file(params.rn_scale_slope))
        } else {
            slope_channel = Channel.value(file("NO_SLOPE"))
        }
        
        // Optional bias for the sigmoid curve
        if (params.rn_scale_bias != null) {    
            bias_channel = channel.value(file(params.rn_scale_bias))
        } else {
            bias_channel = Channel.value(file("NO_BIAS"))
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
                    blacklist_file,
                    control_file
                ).plate_offset
            } else {
                control_dir = Channel.value(file('NO_CONTROL_DIR'))
            }
            
            // Final scaling
            demultiplex_channelstring = scaling_in.map{row -> row[1]}.collect().map({it.unique().join(" ")})

            // Start running only if all the plate offsets have been calculated
            scaling_channel = calculate_scaling_factors(cellpose_in.last(),
                blacklist_file,
                plates,
                manifest_registration_file,
                control_dir.collect(), // ensures this only runs when control dir is done
                demultiplex_channelstring).scaling_factors.first() // Emit the scale factor file as a value channel
            
            // Old way, emit the file text into the channel
            //    .scaling_factors.map { file -> 
            //        file.text
            //}
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
                (row[9][0] == "none") ? "none" : row[9].collect{it -> row[0] + "=" + (it.toInteger() -1).toString()}.join(" ") // mask channels   
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
            finalize_out = finalize(finalize_in, flatfield_out_string, scaling_channel, slope_channel, bias_channel)[0]
            
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
            cellprofiler_out = finalize_and_cellprofiler(finalize_in, flatfield_out_string, scaling_channel, slope_channel, bias_channel)
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