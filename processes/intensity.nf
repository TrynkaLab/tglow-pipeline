#!/usr/bin/env nextflow



// Measure per-object and per-image intensity features on the unscaled finalized images.
// Runs regardless of whether scaling is enabled - feeds calculate_scaling_factors when
// rn_autoscale is set, and is kept around either way to support a future QC step.
process measure_intensity {
    label 'normal'

    conda params.tg_conda_env
    container params.tg_container

    storeDir "$params.rn_publish_dir/scaling/measurements/$outdir"
    scratch params.rn_scratch

    input:
        tuple val(well), path(images, stageAs: "input_images/*")
        val outdir
    output:
        tuple val(well), path("${well.relpath}/*.parquet"), emit: measurements
    script:
        """
        mkdir -p input/${well.relpath}
        ln -s \$(pwd)/input_images/* input/${well.relpath}/

        measure_intensity_features.py \
        --input input \
        --plate ${well.plate} \
        --well ${well.well} \
        --output ./
        """
    stub:
        """
        mkdir -p ${well.relpath}
        touch ${well.relpath}/object_features.parquet
        touch ${well.relpath}/image_features.parquet
        """
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