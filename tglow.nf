#!/usr/bin/env nextflow



// Load workflows

// Stage
include { stage } from './workflows/stage.nf'
include { prepare_manifest; fetch_raw } from './processes/staging.nf'

// Run pipeline
include { run_pipeline } from './workflows/run_pipeline.nf'
include { convertChannelType; parseManifestFlatfields; readBlacklist} from './lib/utils.nf'
include { deconvolute } from './processes/decon.nf'
include { estimate_flatfield; stage_global_flatfield } from './processes/flatfields.nf'
include { register } from './processes/registration.nf'
include { calculate_scaling_factors; calculate_plate_offsets } from './processes/scaling.nf'
include { cellpose } from './processes/segmentation.nf'
include { cellprofiler; finalize; finalize_and_cellprofiler } from './processes/features.nf'

// Run subcell
include { run_subcell } from './workflows/run_subcell.nf'
include { fetch_model; subcell; index_imagedir} from './processes/subcell.nf'


        
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