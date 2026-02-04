#!/usr/bin/env nextflow

include { parseManifestFlatfields } from '../lib/utils.nf'
include { estimate_flatfield; stage_global_flatfield } from '../processes/flatfields.nf'

import ManifestRecord
import RegistrationRecord

workflow flatfield_estimation {
    
    take:
        manifest
        manifest_registration
        blacklist_file
        rn_manifest
        rn_manifest_registration
        bp_run
        bp_channels
        bp_global_flatfield
        plates
    main:
        //------------------------------------------------------------
        // Run estimate_flatfield
        //------------------------------------------------------------
        // If there is no global overide on flatfield channels, get them from manifest
        def run_flatfield = false
        if (bp_channels == null) {
            // Manually parse the manifest to also grab the first record which is hard to do with channels
            flatfield_channels = parseManifestFlatfields(rn_manifest)
                        
            if (flatfield_channels[0].size() > 0) {run_flatfield = true}
            
            // These have the channel already zero indexed
            flatfield_in = Channel.from(flatfield_channels[0])
            flatfield_in_global = Channel.from(flatfield_channels[1])
                   
                            
            if (rn_manifest_registration) {                
                // Create a flat [<plate> <cycle>] channel
                plate_cycle = manifest_registration
                .map{row -> [row.ref_plate, *row.qry_plates]}
                .flatMap { row -> row.withIndex().collect { item, index ->[item, index] }}
                

                // Create a channel per cycle <ref_plate> <cycle> [<plate1> <plateN>]
                // ref_plate here is not the same as ref_plate in the registration, but just indicates where the global flatfield will be saved 
                 per_cycle = plate_cycle
                .groupTuple(by:1)
                .map{row -> tuple(row[0][0], row[1], row[0])}
                                
                // key, cycle, plate, plate(s), channel, index xml
                flatfield_in_global = per_cycle
                .combine(manifest.map{row -> tuple(row.plate, row)}, by:0)
                .flatMap(row -> row[3].bp_channels.collect{channel -> tuple(row[1] + ":" + channel.toString(), row[1], row[2][0], row[2], channel, row[3].index_xml)}) 
                                
                // These have already been converted to zero indexed
                flatfield_in = flatfield_in
                .combine(plate_cycle, by:0)
                .map{row -> tuple(row[3] + ":" + row[1], row[3], row[0], [row[0]],row[1], row[2])} // key, cycle, plate, plate(s), channel, index xml
                
            } else {
                // TODO: this will likely crash
                flatfield_in_global = flatfield_in_global
                .combine(plates)
                .map{row -> tuple( "0:" + row[1], 0, row[0], row[3].split(" "), row[1], row[2])} // key, cycle, plates, channel, index_xml
                
                flatfield_in = flatfield_in
                .combine(plates)
                .map{row -> tuple("0:" + row[1], 0, row[0], row[3].split(" "), row[1], row[2])} // key, cycle, plates, channel, index_xml
            }
        } else {
            flatfield_in = Channel.from(bp_channels)
        }
        
        if (!bp_run) {
            run_flatfield=false
        }
        
        if (run_flatfield) {
            if (bp_global_flatfield) {
                // Runs only one model per plate
                global_flatfield = estimate_flatfield(flatfield_in_global, file(params.rn_image_dir), blacklist_file).flatfield_out
                
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
                
                flatfield_out = stage_global_flatfield(global_flatfield_in).flatfield_out
            } else {
                flatfield_out = estimate_flatfield(flatfield_in, file(params.rn_image_dir), blacklist_file).flatfield_out
            }
            
            // Concat to plate string format // plate : channel = path
            flatfield_out_string = flatfield_out.map{ row -> (
                row[2] + "_ch" + row[4] + "=" + ((row[5] instanceof List) ? row[5].findAll{!it.fileName.name.startsWith("global_refplate").fileName}[0] : row[5].fileName)
            )
            }.collect().map{it.join(" ")}
            
            // Combine into a single channel with files tuple, and plate=file string
            flatfield_out_final = flatfield_out
            .map{row -> file(row[5])}
            .collect()
            .concat(flatfield_out_string)
            .toList()
            .map({ files, string -> tuple(files, string) })
            
        } else {
            flatfield_out_final = Channel.value( tuple( [file('NO_FLATFIELD')], "NO_FLATFIELD" ) )
            //flatfield_out_string = Channel.value(false)
        }
        
        emit:
            flatfield_out = flatfield_out_final
            //flatfield_out_string = flatfield_out_string
         
    
}