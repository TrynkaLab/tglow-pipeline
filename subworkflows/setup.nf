#!/usr/bin/env nextflow

include { index_imagedir } from '../processes/staging.nf'
include { readBlacklist } from '../lib/utils.nf'
include { checkParamsBase } from '../lib/checks.nf'

workflow setup {
    
    take:
        params
    main:
        
        if (params.rn_skip_checks) {
            log.warn("Skipping parameter checks")
        } else {
            checkParamsBase(params)
        }
        //------------------------------------------------------------
        // Read manifests & input files
        //------------------------------------------------------------
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> 
                new ManifestRecord(
                    plate: row.plate,
                    index_xml: row.index_xml,
                    channels: (row.channels == null) ? "none" : row.channels,
                    ff_channels: (row.ff_channels == null) ? "none" : row.ff_channels,
                    cp_nucl_channel: row.cp_nucl_channel,
                    cp_cell_channel: row.cp_cell_channel,
                    dc_channels: row.dc_channels,
                    dc_psfs: row.dc_psfs,
                    mask_channels: (row.mask_channels == null) ? "none" : row.mask_channels
            )
        }
        
        // Build a value channel for the plates in the manifest
        // '<plate_1> <plate_2> <plate_N>'
        plates = manifest.map{row -> row.plate}.collect().map { it.join(' ') }
        
        //------------------------------------------------------------------------
        // Registration manifest
        if (params.rn_manifest_registration != null) {
            manifest_registration = Channel
            .fromPath(params.rn_manifest_registration)
            .splitCsv(header:true, sep:"\t")
            .map { row -> 
                new RegistrationRecord(
                    ref_plate: row.reference_plate,
                    ref_channel: row.reference_channel,
                    qry_plates: row.query_plates,
                    qry_channels: row.query_channels
                )
            }               
        } 
        
        // Registration manifest, if missing just an empty channel
        if (params.rn_manifest_registration == null) {
            //log.info("No registration provided")
            manifest_registration_file = Channel.value(file('NO_REGISTRATION_MANIFEST'))
            manifest_registration = Channel.empty()
        } else {
            manifest_registration_file = Channel.value(file(params.rn_manifest_registration))
        }
            
        //------------------------------------------------------------------------
        // Blacklist channel, if missing just an empty channel
        if (params.rn_blacklist == null) {
            //log.info("No blacklist provided")
            blacklist_file = Channel.value(file('NO_BLACKLIST'))
        } else {
            blacklist_file = Channel.value(file(params.rn_blacklist))
        }
        
        //------------------------------------------------------------------------
        // Control list channel, if missing just an empty channel
        if (params.sc_control_list == null) {
            //log.info("No controlist provided")
            control_file = Channel.value(file('NO_CONTROL_LIST'))
        } else {
            control_file = Channel.value(file(params.sc_control_list))
        }

        //------------------------------------------------------------
        // Prepare per well input channels
        //------------------------------------------------------------

        // Loop over previously generated manifests assuming stage has been run
        if (params.rn_manifest_well == null) {
            manifests_in = index_imagedir(params.rn_image_dir, file(params.rn_image_dir), manifest.map{row -> row.plate}.unique())
        } else {
            manifests_in = Channel.from(params.rn_manifest_well.split(','))
        }

        // Construct the channel on the well level
        // This was needed without the indexing step file(manifest_path)
        well_channel = manifests_in.flatMap{ manifest_path -> manifest_path.splitCsv(header:["well", "row", "col", "plate"], sep:"\t")}
            .map{ row -> new Well(well: row.well, row: row.row, col: row.col, plate: row.plate) }

        // Filter blacklist. Blacklist read into arrat of <plate>:<well>
        if (params.rn_blacklist != null) {
            blacklist = readBlacklist(params.rn_blacklist)
            well_channel = well_channel.filter{ row -> row.key !in blacklist }
        }

        // Filter to specific wells, usefull for testing
        if (params.rn_wells != null) {
            wells = params.rn_wells.split(",")
            well_channel = well_channel.filter{ row -> row.well in wells }
        }

        // Add the plate for easier combining later
        well_channel = well_channel.map{ row -> tuple(row.plate, row) }

    emit:
        manifest=manifest
        plates=plates
        manifest_registration=manifest_registration
        well_channel=well_channel

        // File channels
        manifest_registration_file=manifest_registration_file
        blacklist_file=blacklist_file
        control_file=control_file

}