#!/usr/bin/env nextflow

import Well
import ManifestRecord

// Deconvolute
process deconvolute {
    //label 'gpu_midmem'
    label params.dc_label
    conda params.tg_conda_env
    storeDir "${params.rn_decon_dir}"
    scratch params.rn_scratch

    input:
        tuple val(well), val(manifest), val(psf_string), path(psfs), path(images, stageAs: "input_images/")
    output:
        tuple val(well), val(manifest), path("${well.plate}/${well.row}/${well.col}/")
    script:
        cmd =
        """
        # Workarround as we cannot use variables from the same tuple in stageAs
        mkdir -p input/${well.plate}/${well.row}
        ln -s \$(pwd)/input_images/* input/${well.plate}/${well.row}

        run_richardson_lucy.py \
        --input input/ \
        --plate ${well.plate} \
        --well ${well.well} \
        --psf ${psf_string} \
        --output ./ \
        --clip_max $params.dc_clip_max \
        --niter $params.dc_niter\
        """
        
        if (params.rn_max_project & !params.rn_hybrid) {
            cmd += " --max_project"
        }
        
        if (params.dc_mode) {
            cmd += " --mode $params.dc_mode"
        }
        
        if (params.dc_regularization) {
            cmd += " --regularization $params.dc_regularization"
        }
        
        if (params.dc_psf_subsample_z) {
            cmd += " --psf_subsample_z $params.dc_psf_subsample_z"
        }
        
        if (params.dc_psf_crop_z) {
            cmd += " --psf_crop_z $params.dc_psf_crop_z"
        }
        
        cmd
    stub:
        """
        mkdir -p "${well.relpath}"
        cd "${well.relpath}"
        touch 1.ome.tiff
        """
}
