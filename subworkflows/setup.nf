#!/usr/bin/env nextflow


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
        plates = manifest.map{row -> row[0]}.collect().map { it.join(' ') }
        
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