#!/usr/bin/env nextflow


include { convertChannelType; parseManifestFlatfields; readBlacklist} from '../lib/utils.nf'
include { fetch_model; subcell} from '../processes/subcell.nf'
include { index_imagedir } from '../processes/staging.nf'


workflow run_subcell {
    
    main:
        
        if( !file("$params.rn_publish_dir/processed_images").exists() ) {
            error "$params.rn_publish_dir/processed_images does not exist. Please run '-entry run_pipeline' with --rn_cache_images true"
        }
            
        // Fetch the model
        model_channel = fetch_model(params.sc_model_ref_channels,
                                    params.sc_model, 
                                    Channel.fromPath("${projectDir}/assets/subcell/models/${params.sc_model_ref_channels}/${params.sc_model}/model_config.yaml", checkIfExists: true),
                                    Channel.fromPath("${projectDir}/assets/subcell/models_urls.yaml", checkIfExists: true))
                        .first()
    
    
        // Create a channel from available plates
        plate_channel = Channel
        .fromPath("$params.rn_publish_dir/processed_images/*", type: 'dir')
        .map{ file -> file.getBaseName() }
                
        // Create manifests if they are missing
        manifests_in = index_imagedir("processed_images", "x", plate_channel)
                
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
            
        subcell(well_channel, model_channel, Channel.value(file("${params.rn_publish_dir}/processed_images")))

    
}