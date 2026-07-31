#!/usr/bin/env nextflow

// Determine scaling factors
process calculate_scaling_factors {
    label 'normal'
    
    conda params.tg_conda_env
    container params.tg_container
    
    storeDir "$params.rn_publish_dir/scaling"
   // publishDir "$params.rn_publish_dir/", mode: 'copy'

    input:
        val x
        path blacklist
        val plates
        path registration_manifest
        path control_dir, stageAs: "control_intensities/*"
        val mask_channels
    output:
        path "scaling_factors.txt", emit: scaling_factors
        path "sigmoid_bias.txt", emit: scaling_bias
        path "sigmoid_slope.txt", emit: scaling_slope
        path "raw_scaling_factors.txt"
        path "intensity_summary.tsv"
        path "channel_index_with_scaling.tsv"
        path "histograms"
    script:
        cmd = 
        """
        calculate_scaling_factors.py \
        --output ./ \
        --plate $plates \
        --q1 $params.rn_autoscale_q1 \
        --q2 $params.rn_autoscale_q2\
        """
        
        // Reads the unscaled finalized images (flatfield/registration/max-projection already
        // baked in) rather than raw/decon input - cheaper to read, and lets scaling factors
        // reflect the same geometry cellprofiler/cellcrops will consume. Only reachable when
        // rn_cache_images=true (enforced in checks.nf), since that's what produces this dir.
        cmd += " --input $params.rn_publish_dir/processed_images/unscaled --scale_max 65535"

        // Add optional blacklist
        if (params.rn_blacklist) {
            cmd += " --blacklist $blacklist"
        }
        
        if (control_dir.fileName.name != "NO_CONTROL_DIR") {
             cmd += " --control_dir control_intensities"
        }
        
        // add optional registration manifest
        if (params.rn_manifest_registration) {
            cmd += " --plate_groups $registration_manifest"
        }
        
        if (params.rn_hybrid & mask_channels != "none") {
            cmd += " --mask_channels $mask_channels"
        }
        
        // TMP dummy variable
        //cmd = "echo test > scaling_factors.txt"   
        cmd
    stub:
        """
        touch scaling_factors.txt
        touch raw_scaling_factors.txt
        touch intensity_summary.tsv
        touch channel_index_with_scaling.tsv
        mkdir histograms
        """
}

