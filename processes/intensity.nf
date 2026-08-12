#!/usr/bin/env nextflow



// Measure per-object and per-image intensity features on the unscaled finalized images.
// Runs regardless of whether scaling is enabled - feeds calculate_scaling_factors when
// sc_autoscale is set, and is kept around either way to support a future QC step.
process measure_intensity {
    label params.sc_label

    conda params.rn_conda_env
    container params.rn_container

    publishDir { "$params.rn_publish_dir/rr__features/measurements/$outdir" }
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
        tuple val(plate), path(features, stageAs: "measurements_in/*")
    output:
        tuple val(plate), path("${plate}/")
    script:
        // Workaround as we cannot use variables from the same tuple in stageAs -
        // stage generically, then symlink into the plate-named dir
        // calculate_scaling_factors expects (Nextflow auto-disambiguates the
        // same-named object/image_features.parquet files coming from different wells).
        """
        mkdir -p ${plate}
        ln -s \$(pwd)/measurements_in/* ${plate}/
        echo "Staged \$(ls ${plate}/ | wc -l) files for plate ${plate}"
        """
    stub:
        """
        mkdir -p ${plate}
        touch ${plate}/stub_object_features.parquet
        touch ${plate}/stub_image_features.parquet
        """
}