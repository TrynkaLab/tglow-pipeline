#!/usr/bin/env nextflow

// Determine scaling factors and sigmoid bias/slope from per-plate object/image intensity
// parquet files (produced by measure_intensity, aggregated per-plate by stage_as_plate).
process calculate_scaling_factors {
    label params.sc_label

    conda params.rn_conda_env
    container params.rn_container

    publishDir "$params.rn_publish_dir/rr__scaling"

    input:
        val x
        path plate_dirs, stageAs: "measurements/*"
        path channel_map
        path controls
    output:
        path "scaling_factors.txt", emit: scaling_factors
        path "sigmoid_bias.txt", emit: scaling_bias
        path "sigmoid_slope.txt", emit: scaling_slope
        path "scaling_index.tsv"
    script:
        cmd =
        """
        calculate_scaling_factors.py \
        --input measurements \
        --output ./ \
        --q1 $params.sc_autoscale_q1 \
        --q2 $params.sc_autoscale_q2 \
        --registration_threshold $params.sc_registration_thresh \
        --registration_feature_pattern $params.sc_registration_pattern \
        """

        if (channel_map.name != "NO_CHANNEL_MAP") {
            cmd += " --channel_map $channel_map"
        }

        if (controls.name != "NO_CONTROLS") {
            cmd += " --controls $controls"
        }

        if (params.sc_skip_sigmoid) {
            cmd += " --skip_sigmoid"
        }

        cmd
    stub:
        """
        touch scaling_factors.txt
        touch sigmoid_bias.txt
        touch sigmoid_slope.txt
        touch scaling_index.tsv
        """
}

