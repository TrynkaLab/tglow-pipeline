#!/usr/bin/env nextflow

import ManifestRecord
import RegistrationRecord

workflow setup {
    
    take:
        params
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
                    bp_channels: (row.bp_channels == null) ? "none" : row.bp_channels,
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
        if (params.rn_control_list == null) {
            //log.info("No controlist provided")
            control_file = Channel.value(file('NO_CONTROL_LIST'))
        } else {
            control_file = Channel.value(file(params.rn_control_list))
        }
        
    emit:
        manifest=manifest
        plates=plates
        manifest_registration=manifest_registration
        
        // File channels
        manifest_registration_file=manifest_registration_file
        blacklist_file=blacklist_file
        control_file=control_file
    
}