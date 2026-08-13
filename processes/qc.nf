#!/usr/bin/env nextflow

// Build one "before decon" / "after decon" side-by-side PNG per channel for a single
// sampled well. Reuses stage_cellprofiler.py (the same tool `finalize` uses to apply
// registration and merge cycle channels) against the raw image dir and the decon dir in
// turn, so the comparison is made at the same final, registration-merged channel layout
// the rest of the pipeline uses - without reimplementing that merge logic here.
process qc_decon_samples {
    label params.sc_label

    conda params.rn_conda_env
    container params.rn_container

    scratch params.rn_scratch

    input:
        tuple val(well_id), val(ref_plate), val(qry_plates)
        path raw_image_dir
        path decon_dir
        path registration_dir
        val decon_barrier
    output:
        path "before_after/*.png", emit: images
    script:
        merge_flag = qry_plates ? "--plate_merge " + qry_plates.join(" ") : ""
        registration_flag = (registration_dir.name != "NO_REGISTRATION") ? "--registration_dir ${registration_dir}" : ""
        """
        mkdir -p before after before_after

        stage_cellprofiler.py \
        --input ${raw_image_dir} \
        --output ./before \
        --output_format OME_TIFF \
        --well ${well_id} \
        --plate ${ref_plate} \
        --fields 1 \
        ${merge_flag} \
        ${registration_flag}

        stage_cellprofiler.py \
        --input ${decon_dir} \
        --output ./after \
        --output_format OME_TIFF \
        --well ${well_id} \
        --plate ${ref_plate} \
        --fields 1 \
        ${merge_flag} \
        ${registration_flag}

        make_decon_before_after.py \
        --before ./before \
        --after ./after \
        --plate ${ref_plate} \
        --well ${well_id} \
        --output ./before_after
        """
    stub:
        """
        mkdir -p before_after
        touch before_after/${ref_plate}_${well_id}_ch0_before_after.png
        """
}


// Renders the final self-contained qc_report.html from whichever tabs are available.
// Missing optional inputs are passed as NO_X placeholders (same convention used
// elsewhere in the pipeline) - render_qc_report.py treats a nonexistent path as "tab
// not available" rather than erroring.
process qc_render_report {
    label "tiny"

    conda params.rn_conda_env
    container params.rn_container

    publishDir "$params.rn_publish_dir/rr__qc", mode: 'copy'

    input:
        path plate_dirs, stageAs: "measurements/*"
        path manifest_file
        path blacklist_file
        path registration_manifest_file
        val registration_barrier
        path registration_images_dir
        val show_flatfield
        val flatfield_barrier
        path flatfields_dir
        val ff_global_flatfield
        val ff_params_json
        val show_decon
        path decon_pngs, stageAs: "decon_samples/*"
        val show_scaling
        path scaling_index
    output:
        path "qc_report.html"
    script:
        decon_flag = show_decon ? "--show_decon" : ""
        scaling_flag = show_scaling ? "--show_scaling" : ""
        """
        render_qc_report.py \
        --measurements_dir measurements \
        --manifest ${manifest_file} \
        --output qc_report.html \
        --blacklist ${blacklist_file} \
        --registration_manifest ${registration_manifest_file} \
        --registration_images_dir ${registration_images_dir} \
        --qc_regcor ${params.qc_regcor} \
        --qc_registration_pattern ${params.qc_registration_pattern} \
        --qc_plate_format ${params.qc_plate_format} \
        --qc_n_sample_registration ${params.qc_n_sample_registration} \
        ${show_flatfield ? "--show_flatfield" : ""} \
        --flatfields_dir ${flatfields_dir} \
        ${ff_global_flatfield ? "--ff_global_flatfield" : ""} \
        --ff_params '${ff_params_json}' \
        ${decon_flag} \
        --decon_samples_dir decon_samples \
        ${scaling_flag} \
        --scaling_index ${scaling_index}
        """
    stub:
        """
        touch qc_report.html
        """
}
