#!/usr/bin/env nextflow

process finalize {
    label params.cpr_label
    conda params.cpr_conda_env
    storeDir "$params.rn_publish_dir/processed_images"
    scratch params.rn_scratch

    input:
        tuple val(key), val(well), val(manifest), path(cell_masks, stageAs: "cell_masks/*"), path(nucl_masks, stageAs: "nucl_masks/*"), val(decon_plates), val(merge_plates), path(registration, stageAs:"registration_tmp/*")
        path images, stageAs: "input_images"
        tuple path(basicpy_files), val(basicpy_string)
        path scaling_file
        path slope_file
        path bias_file
    output:
        tuple val(well), path("${well.relpath}/*.tiff")
        path "${well.plate}/channel_indices.tsv"
    script:
    
        // Stage the masks so cellprofiler can access them
        cmd =
        """
        mkdir -p ./masks/${well.relpath}"
        ln -s cell_masks/*  ./masks/${well.relpath}/
        """
        
        if (!nucl_masks[0].name.startsWith("NO_NUCL_MASK")) {
            cmd += "ln -s nucl_masks/* ./masks/${well.relpath}/"
        }
        
        // Stage the registration 
        cmd += 
        """
        mkdir -p registration/${well.plate}/${well.row}/
        ln -s registration_tmp/* registration/${well.plate}/${well.row}/
        """
        
        // Outputs the proccessed images
        cmd += 
        """
        # Stage files
        stage_cellprofiler.py \
        --input input_images \
        --output ./ \
        --output_format OME_TIFF \
        --well ${well.well} \
        --plate ${well.plate} \
        """
        
        if (merge_plates) {
            cmd += " --plate_merge " + merge_plates
        }
        
        if (registration.fileName.name != "NO_REGISTRATION") {
            cmd += " --registration_dir ./registration"
        }
                
        if (basicpy_string) {
            cmd += " --flatfields $basicpy_string"
        }
        
        if ((params.rn_manualscale != null | params.rn_autoscale) & scaling_file.name != "NO_SCALE")  {
            cmd += " --scaling_factors $scaling_file"
        }
        
        if (slope_file.name != "NO_SLOPE")  {
            cmd += " --scaling_slope $slope_file"
        }  
                
        if (bias_file.name != "NO_BIAS")  {
            cmd += " --scaling_bias $bias_file"
        }
        
        if (params.rn_max_project | params.rn_hybrid) {
            cmd += " --max_project --no_zstack"
        }
        
        if (params.rn_hybrid & manifest.mask_channels != "none") {
            cmd += " --mask_dir ./masks"
            cmd += " --mask_pattern *_nucl_mask_*_cp_masks.tiff"
            cmd += " --mask_channels ${manifest.mask_channels.join(' ')}"
        }
        
        if (params.rn_hybrid) {
            cmd += 
            """
            max_project.py \
            --input ./masks \
            --output ./ \
            --well ${well.well} \
            --plate ${well.plate} \
            --pattern *_cell_mask_*_cp_masks.tiff \
            --suffix _cell_mask_d00_ch0_cp_masks.tiff
            """
            
            if (!nucl_masks[0].fileName.name.startsWith("NO_NUCL_MASK")) {
                cmd += 
                """
                max_project.py \
                --input ./masks \
                --output ./ \
                --well ${well.well} \
                --plate ${well.plate} \
                --pattern *_nucl_mask_*_cp_masks.tiff \
                --suffix _nucl_mask_d00_ch0_cp_masks.tiff
                """
            }

        } else {
            cmd += "\nmv ./masks/${well.relpath}/* ./${well.relpath}/"
        }
        
        cmd
    stub:
        """
        mkdir -p ${well.relpath}
        cd ${well.relpath}
        touch ${well.plate}_${well.well}_ch0.tiff
        cd ../../
        touch channel_indices.tsv
        """ 
}

