#!/usr/bin/env nextflow
include { validateParameters } from 'plugin/nf-schema'
include { prepare_manifest; parse_fieldmatrix; fetch_raw } from '../processes/staging.nf'

// Workflow to stage the data from NFS to lustre
workflow stage {
    main:
    
        // Check parameters
        if (params.rn_skip_checks) {
            log.warn("Skipping parameter checks")
        } else {
            // Validate parameters plugin
            validateParameters()
                                
            //------------------------------------------------------------
            if (params.rn_publish_dir == null) {
                error "rn_publish_dir file parameter is required: --rn_publish_dir"
            }
            
            if (params.rn_manifest == null) {
                error "rn_manifest file parameter is required: --rn_manifest"
            }
        }


        // Read manifest
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(row.plate, row.index_xml)}
                
        // Read in the manifest
        if (params.rn_manifest_well == null) {
            
            // Prepare the PE manifest based on the index xml
            prepped_manifest = prepare_manifest(manifest)
            
            // Extract the field matrix from the index json
            parse_fieldmatrix(prepped_manifest.out.metadata)
            
            manifests_in = prepped_manifest.out.manifest
            
        } else {
            manifests_in = Channel.from(params.rn_manifest_well)
        }
        
        well_channel = manifests_in.flatMap{ manifest_path -> file(manifest_path).splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }

        // tuple val(key), val(plate), val(well), val(row), val(col), val(index_xml)
        well_channel = well_channel.map{ row -> tuple(
            row.plate + ":" + row.well,
            row.plate,
            row.well,
            row.row,
            row.col,
            row.index_xml
        ) }
        

        // Filter blacklist
        if (params.rn_blacklist != null) {

            // Read blacklist as list of plate:well
            blacklist=[]
            
            new File(params.rn_blacklist).splitEachLine("\t") {fields ->
                blacklist.add(fields[0] + ":" + fields[1])
            }
            
            log.info("Blacklist consists of items: " + blacklist)
            
            well_channel=well_channel.filter{ row -> row[0] !in blacklist }
            
        }
                
        fetch_raw(well_channel)
        
}