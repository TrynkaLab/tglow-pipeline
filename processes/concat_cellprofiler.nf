#!/usr/bin/env nextflow

// Assigns a short, deterministic index (P1, P2, ...) to every plate in the manifest,
// sorted by plate name so it stays stable regardless of Nextflow's actual (parallel,
// non-deterministic) task-completion order. Used to build globally unique cell/image
// ids in concat_cellprofiler, since CellProfiler's own ImageNumber/ObjectNumber only
// need to be unique within one well's output.
process build_plate_index {
    label "tiny"

    conda params.rn_conda_env
    container params.rn_container

    input:
        path manifest
    output:
        path "plate_index.tsv"
    script:
        "build_plate_index.py --manifest $manifest --output plate_index.tsv"
    stub:
        // A header-only file here would make every downstream .collect() silently emit
        // nothing instead of an empty list (same class of bug already worked around in
        // finalize's stub - see its "Fabricate a plausible channel_indices.tsv" comment),
        // so derive real plate names from the staged manifest instead of a hardcoded guess.
        """
        printf 'plate\\tplate_id\\n' > plate_index.tsv
        tail -n +2 $manifest | cut -f1 | sort -u | awk '{print \$0"\\tP"NR}' >> plate_index.tsv
        """
}

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
