#!/usr/bin/env nextflow

// Determine scaling factors and sigmoid bias/slope from per-plate object/image intensity
// parquet files (produced by measure_intensity, aggregated per-plate by stage_as_plate).
process calculate_scaling_factors {
    label 'normal'

    conda params.tg_conda_env
    container params.tg_container

    storeDir "$params.rn_publish_dir/scaling"

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
        """
        calculate_scaling_factors.py \
        --input measurements \
        --channel_map $channel_map \
        --controls $controls \
        --output ./ \
        --q2 $params.rn_autoscale_q2
        """
    stub:
        """
        touch scaling_factors.txt
        touch sigmoid_bias.txt
        touch sigmoid_slope.txt
        touch scaling_index.tsv
        """
}

