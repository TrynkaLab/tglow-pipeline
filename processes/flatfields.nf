#!/usr/bin/env nextflow

// Estimate_flatfield
process estimate_flatfield {
    label params.ff_label
    
    conda params.rn_conda_env
    container params.rn_container

    storeDir { "${params.rn_publish_dir}/flatfields/${plate}" } //, mode: 'copy'

    input:
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel), val(pe_index)
        path images, stageAs: 'input_images'
        path blacklist
    output:
        path "*${plate}_ch${img_channel}", emit: basicpy_file_out
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel), path("*${plate}_ch${img_channel}"), emit: flatfield_out       
    script:
        cmd =
        """
        run_flatfield_estimation.py \
        --mode $params.ff_mode \
        --input $images \
        --output ./ \
        --nimg $params.ff_nimg \
        --plot \
        --pe_index '$pe_index' \
        --onemodel \
        --channel $img_channel\
        """
        
        if (params.ff_global_flatfield) {
            cmd += " --output_prefix global_refplate_$plate"
            cmd += " --plate " + plates.join(' ')
        } else {
            cmd += " --output_prefix $plate"
            cmd += " --plate $plate"
        }
        
        if (params.rn_max_project) {
            cmd += " --max_project"
        }
        
        if (params.ff_nimg_test) {
            cmd += " --nimg_test $params.ff_nimg_test"
        }
        
        if (params.ff_no_tune) {
            cmd += " --no_tune"
        }   
        
        if (params.ff_merge_n) {
            cmd += " --merge_n $params.ff_merge_n"
        }
        
        if (params.ff_degree) {
            cmd += " --degree $params.ff_degree"
        }
        
        if (params.ff_use_ridge) {
            cmd += " --ridge"
        }
         
        if (params.ff_pseudoreplicates) {
            cmd += " --pseudoreplicates $params.ff_pseudoreplicates"
        }
        
        if (params.ff_pseudoreplicates_test) {
            cmd += " --pseudoreplicates_test $params.ff_pseudoreplicates_test"
        }
        
        if (params.ff_all_planes) {
            cmd += " --all_planes"
        }   
        
        if (params.ff_threshold) {
            cmd += " --threshold"
        }
        
        if (params.ff_autosegment) {
            cmd += " --autosegment"
        }
        
        if (params.rn_blacklist) {
            cmd += " --blacklist $blacklist"
        }
          
        cmd
    stub:
        """
        mkdir ${plate}_ch${img_channel}
        cd ${plate}_ch${img_channel}
        touch flatfield.npy
        """
}

// Copy the global flatfield to per plate folders
process stage_global_flatfield {
    label "tiny"
    
    conda params.rn_conda_env
    container params.rn_container
    
    storeDir { "${params.rn_publish_dir}/flatfields/${plate}/" } //, mode: 'copy'
    //publishDir "${params.rn_publish_dir}/flatfields/${plate}", mode: 'copy'
    
    input:
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel),  path(refdir)
    output:
        path "${plate}_ch${img_channel}", emit: basicpy_file_out
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel),  path("${plate}_ch${img_channel}"), emit: flatfield_out      
    script:
        """
        cp -r $refdir ${plate}_ch${img_channel}
        """
    stub:
        """
        mkdir ${plate}_ch${img_channel}
        cd ${plate}_ch${img_channel}
        touch flatfield.npy
        """
}
