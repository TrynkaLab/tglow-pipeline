#!/usr/bin/env nextflow



process measure_intensity {

    label params.cp_label


    conda params.tg_conda_env
    container params.tg_container

    storeDir "${params.rn_publish_dir}/features/unscaled_intensity"
    scratch params.rn_scratch

    input:
        tuple val(key), val(well),  path(cell_masks, stageAs: "cell_masks/*"), path(nucl_masks, stageAs: "nucl_masks/*"), val(decon_plates), val(merge_plates), path(registration, stageAs:"registration_tmp/*")

    output:
        path val(key), val(well), "features/${well.relpath}/*"
    script:
    cmd = 
    """
    
    """
    cmd
}


// Process takes output of measure_intensity and stages it as a plate-level directory 
// for downstream processing. 
// This is to reduce the length of the filelist NF needs to stage at any given time
process stage_as_plate {
    
    label "tiny"
    
    input:
        tuple val(plate), path(features, stageAs: "${plate}/*")
    output:
        tuple val(plate), path("${plate}/")
    
}