#!/usr/bin/env nextflow

// Run a cellprofiler run
// regular queue
process cellprofiler {
    label params.cpr_label
    conda params.cpr_conda_env
    storeDir "$params.rn_publish_dir/cellprofiler"
    scratch params.rn_scratch

    input:
        tuple val(plate), val(key), val(well), val(row), val(col), path(cell_masks), path(nucl_masks), val(merge_plates), path(registration, stageAs:"registration/*"), val(mask_channels)
        val basicpy_string
        val scaling_string
    output:
        path "features/$plate/$row/$col/*"
    script:
    
        // Stage the masks so cellprofiler can access them
        cmd = "\nmkdir -p ./masks/$plate/$row/$col"
        cmd += "\nmv " + cell_masks.join(" ") + " ./masks/$plate/$row/$col/"
        if (!nucl_masks[0].name.startsWith("NO_NUCL_MASK")) {
            cmd += "\nmv " + nucl_masks.join(" ") + " ./masks/$plate/$row/$col/"
        }
    
        // Outputs the cp files into ./images
        cmd += 
        """
        # Stage files
        stage_cellprofiler.py \
        --output ./images \
        --well $well \
        --plate $plate\
        """
        
        if (params.dc_run) {
            cmd += " --input $params.rn_decon_dir"
        } else {
            cmd += " --input $params.rn_image_dir"
        }
        
        if (merge_plates) {
            cmd += " --plate_merge " + merge_plates
        }
        
        if (registration.fileName.name != "NO_REGISTRATION") {
            cmd += " --registration_dir ./registration"
        }
                
        if (basicpy_string) {
            cmd += " --basicpy_model $basicpy_string"
        }
        
        if ((params.rn_manualscale | params.rn_autoscale) & scaling_string != "none")  {
            cmd += " --scaling_factors $scaling_string"
        }
        
        if (params.rn_max_project | params.rn_hybrid) {
            cmd += " --max_project --no_zstack"
        }
        
        if (params.rn_hybrid & mask_channels != "none") {
            cmd += " --mask_dir ./masks"
            cmd += " --mask_pattern *_nucl_mask_*_cp_masks.tiff"
            cmd += " --mask_channels $mask_channels"
        }
        
        if (params.rn_hybrid) {
            cmd += 
            """
            max_project.py \
            --input ./masks \
            --output ./images \
            --well $well \
            --plate $plate \
            --pattern *_cell_mask_*_cp_masks.tiff \
            --suffix _cell_mask_d00_ch0_cp_masks.tiff
            """
            
            if (!nucl_masks[0].fileName.name.startsWith("NO_NUCL_MASK")) {
                cmd += 
                """
                max_project.py \
                --input ./masks \
                --output ./images \
                --well $well \
                --plate $plate \
                --pattern *_nucl_mask_*_cp_masks.tiff \
                --suffix _nucl_mask_d00_ch0_cp_masks.tiff
                """
            }

        } else {
            cmd += "\nmv ./masks/$plate/$row/$col/* ./images/$plate/$row/$col/"
        }

        // Run cell profiler
        cmd +=
        """
        # Run cellprofiler
        cellprofiler \
        -c \
        -r \
        -o ./features/$plate/$row/$col \
        -i ./images/$plate/$row/$col\
        """
        
        if (params.cpr_plugins) {
            cmd += " --plugins-directory $params.cpr_plugins"
        }
    
        if (params.rn_max_project | params.rn_hybrid) {
           cmd += " -p $params.cpr_pipeline_2d"
        } else {
            cmd += " -p $params.cpr_pipeline_3d"
        }      
        
        // Zip output to save of lustre filecount
        if (!params.cpr_no_zip) {
            cmd +=
            """

            zip -r ./features/$plate/$row/$col/${plate}_${well}.zip ./features/$plate/$row/$col/*.txt
            
            # Cleanup so only zip is staged
            rm ./features/$plate/$row/$col/*.txt
            
            """
        }

        cmd
}