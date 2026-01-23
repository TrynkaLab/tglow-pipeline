#!/usr/bin/env nextflow

// Register
process register {
    //label 'normal'
    label params.rg_label
    
    conda params.tg_conda_env
    container params.tg_container
    
    storeDir "${params.rn_publish_dir}/registration/"
    input:
        tuple val(well), val(registration), path(image_dir, stageAs: "input_images")
    output:
        tuple val(well), val(registration), path("${well.plate}/${well.row}/${well.col}")
    script:
        cmd =
        """
        run_registration.py \
        --input input_images/ \
        --output ./ \
        --well ${well.well} \
        --plate ${well.plate} \
        --plate_merge ${registration.qry_plates.join(" ")} \
        --ref_channel ${registration.ref_channel} \
        --qry_channel ${registration.qry_channels.join(" ")} \
        """
        
        if (params.rg_plot) {
            cmd += " --plot"
        }
        
        if (params.rg_offset_x) {
            cmd += " --offset_x $params.rg_offset_x"
        }
        
        if (params.rg_offset_y) {
            cmd += " --offset_y $params.rg_offset_y"
        }
        
        if (params.rg_mode) {
            cmd += " --mode $params.rg_mode"
        }
        
        if (params.rg_eval) {
            cmd +=
            """ \
            --eval_merge \
            --ref_channel_eval ${registration.ref_channel} \
            --qry_channel_eval ${registration.qry_channels.join(" ")} \
            """
        }
        
        cmd
    stub:
        """
        mkdir -p "${well.relpath}"
        cd "${well.relpath}"
        touch registration.npy
        """
}

