#!/usr/bin/env nextflow


// Deconvolute
process deconvolute {
    //label 'gpu_midmem'
    label params.dc_label
    conda params.tg_conda_env
    storeDir "${params.rn_decon_dir}"
    scratch params.rn_scratch

    input:
        tuple val(plate), val(well), val(row), val(col), val(nucl_channel), val(cell_channel), val(psf_string), path(psfs)
    output:
        tuple val(plate), val(well), val(row), val(col), val(nucl_channel), val(cell_channel), path("$plate/$row/$col")
    script:
        cmd =
        """
        run_richardson_lucy.py \
        --input $params.rn_image_dir \
        --plate $plate \
        --well $well \
        --psf $psf_string \
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
        mkdir -p "$plate/$row/$col"
        cd "$plate/$row/$col"
        touch 1.ome.tiff
        """
}
