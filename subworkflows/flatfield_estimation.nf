#!/usr/bin/env nextflow

include { parseManifestFlatfields } from '../lib/utils.nf'
include { estimate_flatfield; stage_global_flatfield } from '../processes/flatfields.nf'

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
                //log.info("Estimating one global flatfield per cycle")
                
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
                global_flatfield = estimate_flatfield(flatfield_in_global, blacklist_file).flatfield_out
                
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
                flatfield_out = estimate_flatfield(flatfield_in, blacklist_file).flatfield_out
            }
            
            // Concat to plate string format // plate : channel = path
            flatfield_out_string = flatfield_out.map{ row -> (
                row[2] + "_ch" + row[4] + "=" + ((row[5] instanceof List) ? row[5].findAll{!it.fileName.name.startsWith("global_refplate")}[0] : row[5])
            )
            }.collect().map{it.join(" ")}
            
        } else {
            flatfield_out = Channel.empty()
            flatfield_out_string = Channel.value(false)
        }
        
        emit:
            flatfield_out = flatfield_out
            flatfield_out_string = flatfield_out_string
         
    
}