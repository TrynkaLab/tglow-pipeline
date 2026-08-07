#!/usr/bin/env nextflow

process finalize {
    label params.cpr_label
    
    conda params.cpr_conda_env
    container params.rn_container
    
    publishDir "$params.rn_publish_dir/rr__processed_images"
    scratch params.rn_scratch

    input:
        tuple val(key), val(well), val(manifest), path(cell_masks, stageAs: "cell_masks/*"), path(nucl_masks, stageAs: "nucl_masks/*"), val(decon_plates), val(merge_plates), path(registration, stageAs:"registration_tmp/*")
        path images, stageAs: "input_images"
        tuple path(basicpy_files), val(basicpy_string)
        path scaling_file
        path slope_file
        path bias_file
    output:
        tuple val(well), path("{scaled,unscaled}/${well.relpath}/*.tiff"), emit: processed_output
        tuple val(well), path("masks/${well.relpath}/*.tiff"), emit: mask_output
        path "{scaled,unscaled}/${well.plate}/channel_indices.tsv", emit: channel_indices
    script:
        // Whether scaling is actually being applied in this task determines which
        // output subdir the images land in - lets callers get direct-to-scaled output
        // for IO-heavy runs (pass real scaling files) or defer to `rescale` later
        // (pass NO_SCALE/NO_SLOPE/NO_BIAS placeholders).
        apply_scaling = (params.sc_manualscale != null || params.sc_autoscale) && scaling_file.name != "NO_SCALE"
        img_subdir = apply_scaling ? "scaled" : "unscaled"

        // Stage the masks so cellprofiler can access them
        cmd =
        """
        mkdir -p ./masks_stage/${well.relpath}
        ln -s \$(pwd)/cell_masks/*  ./masks_stage/${well.relpath}/
        """

        if (!nucl_masks[0].name.startsWith("NO_NUCL_MASK")) {
            cmd += "ln -s \$(pwd)/nucl_masks/* ./masks_stage/${well.relpath}/"
        }
        
        // Stage the registration 
        cmd += 
        """
        mkdir -p registration/${well.plate}/${well.row}/
        ln -s \$(pwd)/registration_tmp/* registration/${well.plate}/${well.row}/
        """
        
        // Outputs the proccessed images
        cmd += 
        """
        # Stage files
        stage_cellprofiler.py \
        --input input_images \
        --output ./${img_subdir} \
        --output_format OME_TIFF \
        --well ${well.well} \
        --plate ${well.plate} \
        """
        
        if (merge_plates) {
            cmd += " --plate_merge " + merge_plates.join(" ")
        }
        
        if (registration.fileName.name != "NO_REGISTRATION") {
            cmd += " --registration_dir ./registration"
        }
                
        if (basicpy_string && basicpy_string != "NO_FLATFIELD") {
            cmd += " --flatfields $basicpy_string"
        }
        
        if (apply_scaling)  {
            cmd += " --scaling_factors $scaling_file"
        }
        
        if (slope_file.name != "NO_SLOPE")  {
            cmd += " --scaling_slope $slope_file"
        }  
                
        if (bias_file.name != "NO_BIAS")  {
            cmd += " --scaling_bias $bias_file"
        }
        
        if (params.rn_max_project || params.rn_hybrid) {
            cmd += " --max_project --no_zstack"
        }

        if (params.rn_hybrid && manifest.mask_channels != "none") {
            cmd += " --mask_dir ./masks_stage"
            cmd += " --mask_pattern *_nucl_mask_*_cp_masks.tiff"
            cmd += " --mask_channels ${manifest.mask_channels.collect{ well.plate + "=" + it }.join(' ')}"
        }

        if (params.rn_hybrid) {
            cmd +=
            """
            max_project.py \
            --input ./masks_stage \
            --output ./masks \
            --well ${well.well} \
            --plate ${well.plate} \
            --pattern *_cell_mask_*_cp_masks.tiff \
            --suffix _cell_mask_d00_ch0_cp_masks.tiff
            """

            if (!nucl_masks[0].fileName.name.startsWith("NO_NUCL_MASK")) {
                cmd +=
                """
                max_project.py \
                --input ./masks_stage \
                --output ./masks \
                --well ${well.well} \
                --plate ${well.plate} \
                --pattern *_nucl_mask_*_cp_masks.tiff \
                --suffix _nucl_mask_d00_ch0_cp_masks.tiff
                """
            }

        } else {
            cmd +=
            """
            mkdir -p ./masks/${well.relpath}
            rsync --copy-links ./masks_stage/${well.relpath}/* ./masks/${well.relpath}/
            """
        }
        
        cmd
    stub:
        apply_scaling = (params.sc_manualscale != null || params.sc_autoscale) && scaling_file.name != "NO_SCALE"
        img_subdir = apply_scaling ? "scaled" : "unscaled"

        // Fabricate a plausible channel_indices.tsv (real format: see
        // ProcessedImageProvider.build_channel_index in tglow-core) instead of
        // touching an empty file - when registration is enabled, finalize_images'
        // cellcrop_in construction parses this to map original -> registered
        // channel indices, and an empty file makes that combine() always empty,
        // silently dropping every well from cellcrops.
        // Assumes merge plates share the ref plate's channel count, since the
        // real count comes from each plate's actual image dims, which the stub
        // has no image to read.
        ref_channels = (manifest.channels instanceof List) ? manifest.channels : []
        channel_rows = []
        channel_id = 0
        ref_channels.each { ch ->
            channel_rows << [well.plate, well.plate, 1, channel_id, "ch${ch}", ch, "ch${ch}"].join("\t")
            channel_id = channel_id + 1
        }
        if (merge_plates) {
            merge_plates.eachWithIndex { mp, cycleIdx ->
                ref_channels.each { ch ->
                    channel_rows << [well.plate, mp, cycleIdx + 2, channel_id, "ch${channel_id} - ch${ch}", ch, "ch${ch}"].join("\t")
                    channel_id = channel_id + 1
                }
            }
        }
        channel_indices_tsv = (["ref_plate\tplate\tcycle\tchannel\tname\torig_channel\torig_name"] + channel_rows).join("\n")

        """
        mkdir -p ${img_subdir}/${well.relpath}
        touch ${img_subdir}/${well.relpath}/${well.plate}_${well.well}_ch0.tiff
        mkdir -p masks/${well.relpath}
        touch masks/${well.relpath}/${well.plate}_${well.well}_cell_mask_d00_ch0_cp_masks.tiff
        mkdir -p ${img_subdir}/${well.plate}
        printf '%s\\n' "${channel_indices_tsv}" > ${img_subdir}/${well.plate}/channel_indices.tsv
        """
}

// Rescales already-finalized unscaled images (flatfield/registration/max-projection already
// baked in). No raw/decon/multi-plane reads happen here, so this is cheap relative to `finalize`.
process rescale {
    label params.cpr_label

    conda params.rn_conda_env
    container params.rn_container

    publishDir "$params.rn_publish_dir/rr__processed_images"
    scratch params.rn_scratch

    input:
        tuple val(well), path(images, stageAs: "input_images/*")
        path scaling_file
        path slope_file
        path bias_file
    output:
        tuple val(well), path("scaled/${well.relpath}/*.ome.tiff"), emit: processed_output
    script:
        cmd =
        """
        # Workaround as we cannot use variables from the same tuple in stageAs
        mkdir -p unscaled/${well.relpath}
        ln -s \$(pwd)/input_images/* unscaled/${well.relpath}/

        rescale_images.py \
        --input ./unscaled \
        --output ./scaled \
        --plate ${well.plate} \
        --well ${well.well} \
        --scaling_factors $scaling_file \
        """

        if (slope_file.name != "NO_SLOPE") {
            cmd += " --scaling_slope $slope_file"
        }

        if (bias_file.name != "NO_BIAS") {
            cmd += " --scaling_bias $bias_file"
        }

        cmd
    stub:
        """
        mkdir -p scaled/${well.relpath}
        touch scaled/${well.relpath}/${well.plate}_${well.well}_ch0.ome.tiff
        """
}

process cellcrops {
    label params.cpr_label
    
    conda params.rn_conda_env
    container params.rn_container

    publishDir "$params.rn_publish_dir/rr__cellcrops"
    scratch params.rn_scratch

    input:
        tuple val(well), val(registration),  path(images, stageAs: "input_images/"), path(masks, stageAs: "input_masks/")
    output:
        tuple val(well), path("${well.relpath}/*.h5"), emit: h5
        tuple val(well), path("${well.relpath}/*.parquet"), emit: index
    script:

    cmd = """
    mkdir -p input/${well.relpath}
    ln -s \$(pwd)/input_images/* input/${well.relpath}

    mkdir -p input_masks/${well.relpath}
    ln -s \$(pwd)/input_masks/* input_masks/${well.relpath}

    run_cellsampling.py \
    --input input \
    --cell_mask_dir input_masks \
    --nucl_mask_dir input_masks \
    --output ./ \
    --plate ${well.plate} \
    --well ${well.well} \
    --max_per_field $params.rn_max_per_field \
    """
    if (registration != null) {
        cmd += "--ref_channel ${registration.ref_channel} --qry_channels ${registration.qry_channels.join(" ")}"
    }

    cmd
    stub:
        """
        mkdir -p ${well.relpath}
        touch ${well.relpath}/1.h5
        touch ${well.relpath}/1.parquet
        """
}

process index_cellcrops {
    label params.cpr_label

    conda params.rn_conda_env
    container params.rn_container

    publishDir "$params.rn_publish_dir/rr__cellcrops"

    input:
        val previous_completed
        path input, stageAs: "input_cellcrops"
    output:
        path "cellcrop_index.parquet"
    script:
    """
    python3 <<'PYEOF'
import glob

import pandas as pd

paths = sorted(glob.glob("input_cellcrops/**/*.parquet", recursive=True))
df = pd.concat((pd.read_parquet(p) for p in paths), ignore_index=True) if paths else pd.DataFrame()
df.to_parquet("cellcrop_index.parquet")
PYEOF
    """
    stub:
        """
        touch cellcrop_index.parquet
        """
}

// Backstop: aggregates per-field CSVs into a single gzipped CSV, for use only when
// cellcrops is run with --csv_per_field (not exposed as a pipeline param - see
// run_cellsampling.py's --csv_per_field flag).
process index_cellcrops_csv {
    label params.cpr_label

    conda params.rn_conda_env
    container params.rn_container

    publishDir "$params.rn_publish_dir/rr__cellcrops"

    input:
        val previous_completed
        path input, stageAs: "input_cellcrops"
    output:
        path "cellcrop_index.csv.gz"
    script:
    """
    # Get the header from the first CSV file
    find input_cellcrops/ -name "*.csv" -type f -print0 | head -z -n1 | xargs -0 head -n1 > cellcrop_index.csv

    # Append all CSV files without their headers
    find input_cellcrops/ -name "*.csv" -type f -print0 | xargs -0 tail -q -n+2 >> cellcrop_index.csv

    gzip -f cellcrop_index.csv
    """
    stub:
        """
        touch cellcrop_index.csv.gz
        """
}