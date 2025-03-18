
workflow run_subcell {
    
    main:
    
        // Fetch the model
        model_channel = fetch_model("rygb", "mae_contrast_supcon_model")
    
        // Create manifests if they are missing
        manifests_in = index_imagedir("proccessed_images")
        
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
    
        subcell(well_channel, model_channel, file("${params.rn_publish_dir}/proccessed_images"))

    
}