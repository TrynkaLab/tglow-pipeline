#!/usr/bin/env nextflow

// Estimate_flatfield
// normal queue
process estimate_flatfield {
    //label 'himem'
    label params.bp_label
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/flatfields/${plate}" //, mode: 'copy'

    input:
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel), val(pe_index)
        path blacklist
    output:
        path "*${plate}_ch${img_channel}", emit: basicpy_file_out
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel), path("*${plate}_ch${img_channel}"), emit: flatfield_out       
    script:
        cmd =
        """
        run_flatfield_estimation.py \
        --mode $params.bp_mode \
        --input $params.rn_image_dir \
        --output ./ \
        --nimg $params.bp_nimg \
        --plot \
        --pe_index '$pe_index' \
        --onemodel \
        --channel $img_channel\
        """
        
        if (params.bp_global_flatfield) {
            cmd += " --output_prefix global_refplate_$plate"
            cmd += " --plate " + plates.join(' ')
        } else {
            cmd += " --output_prefix $plate"
            cmd += " --plate $plate"
        }
        
        if (params.rn_max_project) {
            cmd += " --max_project"
        }
        
        if (params.bp_nimg_test) {
            cmd += " --nimg_test $params.bp_nimg_test"
        }
        
        if (params.bp_no_tune) {
            cmd += " --no_tune"
        }   
        
        if (params.bp_merge_n) {
            cmd += " --merge_n $params.bp_merge_n"
        }
        
        if (params.bp_degree) {
            cmd += " --degree $params.bp_degree"
        }
        
        if (params.bp_use_ridge) {
            cmd += " --ridge"
        }
         
        if (params.bp_pseudoreplicates) {
            cmd += " --pseudoreplicates $params.bp_pseudoreplicates"
        }
        
        if (params.bp_pseudoreplicates_test) {
            cmd += " --pseudoreplicates_test $params.bp_pseudoreplicates_test"
        }
        
        if (params.bp_all_planes) {
            cmd += " --all_planes"
        }   
        
        if (params.bp_threshold) {
            cmd += " --threshold"
        }
        
        if (params.bp_autosegment) {
            cmd += " --autosegment"
        }
        
        if (params.rn_blacklist) {
            cmd += " --blacklist $blacklist"
        }
          
        cmd
}

// Copy the global flatfield to per plate folders
process stage_global_flatfield {
    label "tiny"
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/flatfields/${plate}/" //, mode: 'copy'

    input:
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel),  path(refdir)
    output:
        path "${plate}_ch${img_channel}", emit: basicpy_file_out
        tuple val(key), val(cycle), val(plate), val(plates), val(img_channel),  path("${plate}_ch${img_channel}"), emit: flatfield_out      
    script:
    """
    cp -r $refdir ${plate}_ch${img_channel}
    """
}
