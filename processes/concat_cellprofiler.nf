#!/usr/bin/env nextflow

// Aggregates every well's CellProfiler zip for a plate into one cell-level and one
// image-level parquet file. Requires cpr_no_zip=false (enforced in checks.nf) - each
// well's zip is already named ${plate}_${well}.zip, so staging many wells' zips into
// one task is collision-free; loose per-well .txt files would collide on filename.
process concat_cellprofiler {
    label params.cpr_label

    conda params.rn_conda_env
    container params.rn_container

    publishDir "$params.rn_publish_dir/rr__features/cellprofiler"
    scratch params.rn_scratch

    input:
        tuple val(plate), val(plate_id), path(zips, stageAs: "zips_in/*")
    output:
        tuple val(plate), path("${plate}_cells.parquet"), path("${plate}_image.parquet")
    script:
        """
        concat_cellprofiler.py \
        --input zips_in \
        --plate ${plate} \
        --plate_id ${plate_id} \
        --merge_strategy ${params.cpr_merge_strategy} \
        --parent_col ${params.cpr_parent_col} \
        --child_pattern '${params.cpr_child_pattern}' \
        --output_cells ${plate}_cells.parquet \
        --output_image ${plate}_image.parquet
        """
    stub:
        """
        touch ${plate}_cells.parquet
        touch ${plate}_image.parquet
        """
}
