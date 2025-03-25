#!/usr/bin/env nextflow

// Register
process register {
    //label 'normal'
    label params.rg_label
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/registration/"
    input:
        tuple val(plate), val(well), val(row), val(col), val(reference_channel), val(query_plates), val(query_channels)
    output:
        tuple val(plate), val(well), val(row), val(col), val(query_plates), path("${plate}/${row}/${col}")
    script:
        cmd =
        """
        run_registration.py \
        --input $params.rn_image_dir \
        --output ./ \
        --well $well \
        --plate $plate \
        --plate_merge $query_plates \
        --ref_channel $reference_channel \
        --qry_channel $query_channels\
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
            --ref_channel_eval $reference_channel \
            --qry_channel_eval $query_channels
            """
        }
        
        cmd

}

