#!/usr/bin/env nextflow

include { convertChannelType } from './lib/utils.nf'

// Prepare a manfiest
process prepare_manifest {
    label 'tiny_img'
    conda params.tg_conda_env
    storeDir "${params.rn_image_dir}/${plate}"
    
    input:
        tuple val(plate), val(index_xml)
    output:
        path "manifest.tsv", emit: manifest
        tuple path("Index.*xml"), path("Index.json"), path("acquisition_info.txt"), emit: metadata
    script:
        cmd =
        """
        parse_xml.py \
        --input_file '$index_xml' \
        --output_path ./ \
        --to_manifest
        
        cp '$index_xml' ./
        """
        
        //def imgdir = new File(index_xml).parent
        //def assaydir = new File(imgdir.text + "/../AssayLayout")
        
        //if (assaydir.exists()) {
        //    cmd += 
        //    """
        //    cp -r '$assaydir' ./
        //    """
        //}
        
        // Only keep the first few lines 
        //if (params.rn_testmode) {
        //    cmd += 
        //    """
        //    head -n 3 manifest.tsv > tmp
        //    mv tmp manifest.tsv
        //    """
        //}
        cmd
    
    //manifest = params.tglow_image_dir + "/" + plate + "/" + "manifest.tsv"
}

// Fetches raw data from NFS and recodes into new OME file structure
// imaging queue
process fetch_raw {
    scratch params.rn_scratch
    label params.st_label
    conda params.tg_conda_env
    storeDir "${params.rn_image_dir}"

    input:
        //tuple val(well), val(row), val(col), val(plate), val(index_xml)
        tuple val(key), val(plate), val(well), val(row), val(col), val(index_xml)

    output:
        path "$plate/$row/$col"
    script:
        """
        convert_pe_raw.py \
        --input_file '$index_xml' \
        --output_path ./ \
        --well $well
        
        md5sum $plate/$row/$col/*.ome.tiff > $plate/$row/$col/CHECKSUMS.txt
        """         
}

// Basicpy
// normal queue
process basicpy {
    //label 'himem'
    label params.bp_label
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/basicpy/${plate}" //, mode: 'copy'

    input:
        tuple val(plate), val(img_channel)
        path blacklist
    output:
        path "${plate}_ch${img_channel}", emit: basicpy_file_out
        tuple val(plate), val(img_channel),  path("${plate}_ch${img_channel}"), emit: basicpy_out       
    script:
        cmd =
        """
        run_basicpy.py \
        --input $params.rn_image_dir \
        --output ./ \
        --output_prefix $plate \
        --plate $plate \
        --nimg $params.bp_nimg \
        --plot \
        --channel $img_channel\
        """
        
        if (params.rn_max_project) {
            cmd += " --max_project"
        }
            
        if (params.bp_no_tune) {
            cmd += " --no_tune"
        }   
        
        if (params.bp_merge_n) {
            cmd += " --merge_n $params.bp_merge_n"
        }
        
        if (params.bp_pseudoreplicates) {
            cmd += " --pseudoreplicates $params.bp_pseudoreplicates"
        }
        
        if (params.bp_all_planes) {
            cmd += " --all_planes"
        }   
        
        if (params.bp_threshold) {
            cmd += " --threshold"
        }
        
        if (params.bp_autosegment) {
            cmd += " --autosegment"
        }
        
        if (params.rn_blacklist) {
            cmd += " --blacklist $blacklist"
        }
          
        cmd
}

// Cellpose
process cellpose {
    label params.cp_label
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/masks/"
    scratch params.rn_scratch

    input:
        tuple val(plate), val(well), val(row), val(col), val(nucl_channel), val(cell_channel)
    output:
        tuple val(plate), val(well), val(row), val(col), path("${plate}/${row}/${col}/*_cell_mask*_ch${cell_channel}*.*"), path("${plate}/${row}/${col}/*_nucl_mask*_ch${nucl_channel}*.*"), emit: cellpose_out //, path("${plate}/${row}/${col}/*_nucl_mask*.tif"), emit: cellpose_out
        
    script:
        cmd =
        """
        run_cellpose.py \
        --output ./ \
        --plate $plate \
        --well $well \
        --cell_channel $cell_channel \
        --gpu \
        --diameter $params.cp_cell_size \
        --model $params.cp_model\
        """
        
        if (params.dc_run) {
            cmd += " --input $params.rn_decon_dir"
        } else {
            cmd += " --input $params.rn_image_dir"
        }
        
        if (nucl_channel >= 0) {
            cmd +=
            """ \
            --nucl_channel $nucl_channel  \
            --diameter_nucl $params.cp_nucl_size\
            """
        }    
        
        if (params.cp_min_cell_area) {
            cmd += " --min_cell_area $params.cp_min_cell_area"
        }
        
        if (params.cp_min_nucl_area) {
            cmd += " --min_nucl_area $params.cp_min_nucl_area"
        }
        
        if (params.cp_cell_flow_thresh) {
            cmd += " --cell_flow_thresh $params.cp_cell_flow_thresh"
        }
        
        if (params.cp_nucl_flow_thresh) {
            cmd += " --nucl_flow_thresh $params.cp_nucl_flow_thresh"
        }
        
        if (params.cp_cell_prob_threshold) {
            cmd += " --cell_prob_threshold $params.cp_cell_prob_threshold"
        }
        
        if (params.cp_nucl_prob_threshold) {
            cmd += " --nucl_prob_threshold $params.cp_nucl_prob_threshold"
        }
        
        if (params.rn_max_project & !params.rn_hybrid) {
            cmd += " --no_3d"
        }
        
        if (params.cp_dont_use_nucl_for_declump) {
            cmd += " --dont_use_nucl_for_declump"
        }
        
        if (params.cp_cell_power) {
            cmd += " --cell_power $params.cp_cell_power"
        }
        
        if (params.cp_nucl_power) {
            cmd += " --nucl_power $params.cp_nucl_power"
        }
        
        if (params.cp_dont_post_process) {
            cmd += " --dont_post_process"
        }
        
        if (params.cp_downsample) {
            cmd += " --downsample $params.cp_downsample"
        }
        
        // Add a fake nucleus channel because nextflow doesn't play nicely with 
        // tuples and mulptiple input files, otherwise making sure wells and channels are 
        // matched turns into a pain
        // If its stupid and it works, it is not stupid, hopefully in future, nextflow
        // will deal better with optional files or just accept null for file objects
        if (nucl_channel < 0) {
            cmd += 
            """
            
            touch ${plate}/${row}/${col}/NO_NUCL_MASK_nucl_mask_ch${nucl_channel}_dummy.tif
            """
        
        }
        cmd

}

// Register
process register {
    //label 'normal'
    label params.rg_label
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/registration/"
    input:
        tuple val(plate), val(well), val(row), val(col), val(reference_channel), val(query_plates), val(query_channels)
    output:
        tuple val(plate), val(well), val(row), val(col), val(query_plates), path("${plate}")
    script:
        cmd =
        """
        run_registration.py \
        --input $params.rn_image_dir \
        --output ./ \
        --well $well \
        --plate $plate \
        --plate_merge $query_plates \
        --ref_channel $reference_channel \
        --qry_channel $query_channels\
        """
        
        if (params.rg_plot) {
            cmd += " --plot"
        }
        
        if (params.rg_eval) {
            cmd +=
            """ \
            --eval_merge \
            --ref_channel_eval $reference_channel \
            --qry_channel_eval $query_channels
            """
        }
        
        cmd

}

// Deconvolute
process deconvolute {
    //label 'gpu_midmem'
    label params.dc_label
    conda params.tg_conda_env
    storeDir "${params.rn_decon_dir}"
    scratch params.rn_scratch

    input:
        tuple val(plate), val(well), val(row), val(col), val(nucl_channel), val(cell_channel), val(psf_string), path(psfs)
    output:
        tuple val(plate), val(well), val(row), val(col), val(nucl_channel), val(cell_channel), path("$plate/$row/$col")
    script:
        cmd =
        """
        run_richardson_lucy.py \
        --input $params.rn_image_dir \
        --plate $plate \
        --well $well \
        --psf $psf_string \
        --output ./ \
        --clip_max $params.dc_clip_max \
        --niter $params.dc_niter\
        """
        
        if (params.rn_max_project & !params.rn_hybrid) {
            cmd += " --max_project"
        }
        
        if (params.dc_mode) {
            cmd += " --mode $params.dc_mode"
        }
        
        if (params.dc_regularization) {
            cmd += " --regularization $params.dc_regularization"
        }
        
        if (params.dc_psf_subsample_z) {
            cmd += " --psf_subsample_z $params.dc_psf_subsample_z"
        }
        
        if (params.dc_psf_crop_z) {
            cmd += " --psf_crop_z $params.dc_psf_crop_z"
        }
        
        cmd
}

// Run a cellprofiler run
// regular queue
process cellprofiler {
    label params.cpr_label
    conda params.cpr_conda_env
    storeDir "$params.rn_publish_dir/cellprofiler"
    scratch params.rn_scratch

    input:
        tuple val(plate), val(key), val(well), val(row), val(col), path(cell_masks), path(nucl_masks), val(merge_plates), path(registration, stageAs:"registration/*"), val(basicpy_string), val(mask_channels)
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
        //else {
            // This was the old way of doing things, now directly staged in mask dir
            //cmd += "\nmv " + cell_masks.join(" ") + " ./images/$plate/$row/$col/"
            //if (!nucl_masks[0].name.startsWith("NO_NUCL_MASK")) {
            //    cmd += "\nmv " + nucl_masks.join(" ") + " ./images/$plate/$row/$col/"
            //}
            // Move raw masks for cellprofiler
        //    cmd += "\nmv ./masks/$plate/$row/$col/* ./images/$plate/$row/$col/"
        //}
         
        //cmd += "\ntouch done.tsv"     
        //cmd

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

// Determine scaling factors
process calculate_scaling_factors {
    label 'normal'
    conda params.tg_conda_env
    storeDir "$params.rn_publish_dir/"
   // publishDir "$params.rn_publish_dir/", mode: 'copy'

    input:
        val x
        path blacklist
        val plates
        path registration_manifest
    output:
        path "scaling_factors.txt", emit: scaling_factors
        path "intensity_summary.tsv"
    script:
        cmd = 
        """
        calculate_scaling_factors.py \
        --plate $plates \
        --q1 $params.rn_autoscale_q1 \
        --q2 $params.rn_autoscale_q2\
        """
        
        if (params.dc_run) {
            // Not sure why I thought the scale_max was needed. The intensity files match the 16 bit already, so this
            // is not needed.
            //cmd += " --input $params.rn_decon_dir --scale_max $params.dc_clip_max"
            cmd += " --input $params.rn_decon_dir --scale_max 65535"
        } else {
            cmd += " --input $params.rn_image_dir --scale_max 65535"
        }
        
        // Add optional blacklist
        if (params.rn_blacklist) {
            cmd += " --blacklist $blacklist"
        }
        
        // add optional registration manifest
        if (params.rn_manifest_registration) {
            cmd += " --plate_groups $registration_manifest"
        }
        
        // TMP dummy variable
        //cmd = "echo test > scaling_factors.txt"   
        cmd
}


// Workflow to stage the data from NFS to lustre
workflow stage {
    main:
        //------------------------------------------------------------
        if (params.rn_publish_dir == null) {
            error "rn_publish_dir file parameter is required: --rn_publish_dir"
        }
        
        if (params.rn_manifest == null) {
            error "rn_manifest file parameter is required: --rn_manifest"
        }

        // Set runtime defaults, these are overidden when specified on commandline
        params.rn_image_dir = params.rn_publish_dir + "/images"
        params.rn_decon_dir = params.rn_publish_dir + "/decon"
    
        // Read manifest
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(row.plate, row.index_xml)}
        
        //manifest.view()
        
        // Read in the manifest
        if (params.rn_manifest_well == null) {
            manifests_in = prepare_manifest(manifest).manifest
        } else {
            manifests_in = Channel.from(params.rn_manifest_well)
        }
        
        //manifests_in.view()
        well_channel = manifests_in.flatMap{ manifest_path -> file(manifest_path).splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }

        // tuple val(key), val(plate), val(well), val(row), val(col), val(index_xml)
        well_channel = well_channel.map(row -> {tuple(
            row.plate + ":" + row.well,
            row.plate,
            row.well,
            row.row,
            row.col,
            row.index_xml
        )})
        

        // Filter blacklist
        if (params.rn_blacklist != null) {

            // Read blacklist as list of plate:well
            blacklist=[]
            
            new File(params.rn_blacklist).splitEachLine("\t") {fields ->
                blacklist.add(fields[0] + ":" + fields[1])
            }
            
            log.info("Blacklist consists of items: " + blacklist)
            
            well_channel=well_channel.filter(row -> {
                row[0] !in blacklist       
            })
            
        }
        
        //well_channel.view()
        
        fetch_raw(well_channel)
        
}

// Main workflow
workflow run_pipeline {

    main:
        //------------------------------------------------------------
        if (params.rn_publish_dir == null) {
            error "rn_publish_dir file parameter is required: --rn_publish_dir"
        }
        
        if (params.rn_manifest == null) {
            error "rn_manifest file parameter is required: --rn_manifest"
        }

        // Set runtime defaults, these are overidden when specified on commandline
        params.rn_image_dir = params.rn_publish_dir + "/images"
        params.rn_decon_dir = params.rn_publish_dir + "/decon"
        
        //------------------------------------------------------------
        // Read manifest
        manifest = Channel.fromPath(params.rn_manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(
            row.plate,
            row.index_xml,
            tuple(row.channels.split(',')), 
            tuple(row.bp_channels.split(',')),
            row.cp_nucl_channel,
            row.cp_cell_channel,
            row.dc_channels,
            row.dc_psfs,
            row.mask_channels,
            row.scaling_factors)}
        //manifest.view()
        
        // Build a value channel for the plates in the manifest
        // '<plate_1> <plate_2> <plate_N>'
        plates_channel = manifest
        .map{row -> row[0]}
        .collect()
        .map { it.join(' ') }


        // Build the value vhannel of manual scaling factors
        // '<plate_1>=<scale_1> <plate_2>=<scale_2> <plate_N>=<scale_N>'
        // TODO: nextflowify this by remaping this from manifest channel
        if (params.rn_manualscale & !params.rn_autoscale) {
            def csvFile = new File(params.rn_manifest)
            def csvData = csvFile.readLines()
            def scaling_string = null

            def plate_scale = []
            for (int i = 1; i < csvData.size(); i++) {
                def curLine = csvData[i].split('\t')
                
                if ((curLine[9] != "none") & (curLine[9] != null)) {
                    def curChannels = curLine[2].split(',')
                    for (channel in curChannels) {
                        curChannel = channel.toInteger()-1
                        plate_scale << curLine[0] + "_ch" + curChannel + "=" + curLine[9].split(",")[curChannel]
                    }
                }
            }   
            
            scaling_string = plate_scale.join(" ")
            log.info("Read scaling factors from manifest: " + scaling_string)
            
            scaling_channel = Channel.value(scaling_string)
        } else {
            scaling_channel = Channel.value("none")
        }
        
        // Blacklist channel, if missing just an empty channel
        if (params.rn_blacklist == null) {
            log.info("No blacklist provided")
            blacklist_channel = Channel.value(file('NO_BLACKLIST'))
        } else {
            blacklist_channel = Channel.value(file(params.rn_blacklist))
        }
    
        // Registration manifest, if missing just an empty channel
        if (params.rn_manifest_registration == null) {
            log.info("No registration provided")
            manifest_registration_channel = Channel.value(file('NO_REGISTRATION_MANIFEST'))
        } else {
            manifest_registration_channel = Channel.value(file(params.rn_manifest_registration))
        }
    
        //------------------------------------------------------------
        // Run basicpy
        
        // If there is no global overide on basicpy channels, get them from manifest
        // TODO: nextflowify this by remaping this from manifest channel
        
        def run_basicpy = false
        if (params.bp_channels == null) {
            def csvFile = new File(params.rn_manifest)
            def csvData = csvFile.readLines()

            def plate_channel = []
            // Skip header
            // Substract one from the channel is python is 0 indexed
            for (int i = 1; i < csvData.size(); i++) {
                def curLine = csvData[i].split('\t')
                
                if (curLine[3] != "none") {
                    def curChannels = curLine[3].split(',')
                    for (channel in curChannels) {
                        plate_channel << tuple(curLine[0], channel.toInteger()-1)
                    }
                }
            }
            
            if (plate_channel.size() > 0) {
                run_basicpy = true
            }
               
            basicpy_in = Channel.from(plate_channel)
        } else {
            basicpy_in = Channel.from(params.bp_channels)
        }
        
        if (run_basicpy) {
            log.info("Running basicpy")
        } else {
            log.info("No manifest entries or --bp_channels provided, so skipping basicpy")
        }
         
        //basicpy_channels = basicpy_channels.combine(manifest, by:0)
        //basicpy_in.view()
        basicpy_out = basicpy(basicpy_in, blacklist_channel).basicpy_out
        
        //------------------------------------------------------------
        // Loop over previously generated manifests assuming stage has been run
        // TODO: nextflowify this by remaping this from manifest channel
        if (params.rn_manifest_well == null) {
            // code that reads paths available manfiests from previous stage
            //manifests_in = Channel.fromPath("${params.rn_image_dir}/*/manifest.tsv")
            
            // Read which plates should be run, in the case the manifest has been edited
            // after running stage to run not all plates
            plates = []
            new File(params.rn_manifest).eachLine{string, index -> 
                if (index > 1) {
                    plates.add(string.split("\t")[0])
                }
            }
            
            manifest_paths = []
            
            for (plate in plates) {
                manifest_paths += "${params.rn_image_dir}/" + plate + "/manifest.tsv"
            }
            
            log.info("Considering plate level manifests: " + manifest_paths)
            manifests_in = Channel.from(manifest_paths)
        } else {
            manifests_in = Channel.from(params.rn_manifest_well.split(','))
        }
        
        // Construct the channel on the well level
        well_channel = manifests_in
            .flatMap{ manifest_path -> file(manifest_path)
            .splitCsv(header:["well", "row", "col", "plate", "index_xml"], sep:"\t") }
                                 
        // Filter blacklist
        if (params.rn_blacklist != null) {
            // Read blacklist as list of plate:well
            blacklist=[]
            
            new File(params.rn_blacklist).splitEachLine("\t") {fields ->
                blacklist.add(fields[0] + ":" + fields[1])
            }
            
            log.info("Blacklist consists of " + blacklist.size() + " items")
            //log.info("Blacklist consists of items: " + blacklist)
            
            well_channel=well_channel.filter(row -> {
                (row.plate + ":" + row.well) !in blacklist       
            })
            
        }
        
        // Filter to specific wells, usefull for testing
        if (params.rn_wells != null) {
            wells = params.rn_wells.split(",")
            log.info("Selecting " + wells.size() + " wells: " + wells)

            well_channel=well_channel.filter(row -> {row.well in wells})
            //well_channel.view()
        }                  
    
        // Re-order the well channel for later merging
        well_in = well_channel.map{ row -> tuple(
            row.plate,
            row.well,
            row.row,
            row.col
        )}
                
        //------------------------------------------------------------
        // Deconvolute
        if (params.dc_run) {
            decon_in = well_in.combine(manifest, by: 0)
            .filter(row -> {row[9] != "none"})
            .map{ row -> tuple(
                row[0], // plate
                row[1], // well
                row[2], // row
                row[3], // col
                convertChannelType(row[7]), // nucl_channel
                convertChannelType(row[8]), // cell_channel
                row[9].split(',').collect{it -> (it.toInteger() -1).toString() + "=" + new File(row[10].split(',')[it.toInteger() -1]).getName()}.join(" "),  // dc_channels
                row[10].split(',').collect(it -> (file(it))) //dc_psfs
            )}      
            
            // Deconvolute
            cellpose_in = deconvolute(decon_in)
            .filter(row -> {row[5] != "none"})
            .map{ row -> tuple(
                row[0], // plate
                row[1], // well
                row[2], // row
                row[3], // col
                row[4], // nucl_channel
                row[5] // cell_channel
            )}            
        } else {
            // Append the things from the manifest needed for cellpose
            cellpose_in = well_in.combine(manifest, by: 0)
            .filter(row -> {row[8] != "none"})
            .map{ row -> tuple(
                row[0], // plate
                row[1], // well
                row[2], // row
                row[3], // col
                convertChannelType(row[7]), // nucl_channel
                convertChannelType(row[8]) // cell_channel
            )}
        }
        
        // Get scaling string
        // When all deconvelution is done, or all data is staged, calculate the scaling factors
        // which map the images onto the full dynamic range of the 16 bit uint
        if (params.rn_autoscale) {
            scaling_channel = calculate_scaling_factors(cellpose_in.last(),
            blacklist_channel,
            plates_channel,
            manifest_registration_channel).scaling_factors.map { file -> 
                file.text
            }
            //.first() this used to be at the end, but removed it as appreantly redundant
            //scaling_channel.view()
        }

        //------------------------------------------------------------
        // Register
        if (params.rn_manifest_registration != null) {
        
            manifest_registration = Channel
            .fromPath(params.rn_manifest_registration)
            .splitCsv(header:true, sep:"\t")
            .map{row -> tuple(
                row.reference_plate,
                row.reference_channel,
                row.query_plates,
                row.query_channels
            )}

            // Append the plate info from the manifest if it exists
            // Use a join so the plates to register are discarded
            // Use -1 as the channels are 0 indexed
            registration_in = well_in
            .combine(manifest_registration, by:0)
            .map{ row -> tuple(
                row[0], // plate
                row[1], // well
                row[2], // row
                row[3], // col
                row[4].toInteger()-1, // reference_channel
                row[5].split(',').join(" "), // query_plates
                row[6].split(',').collect{it -> (it.toInteger() -1).toString()}.join(" ")  // query_channels
            )} 
                 
            // Run registration
            registration_out = register(registration_in)
                          
            // Filter cellpose channel to run reference plates only 
            cellpose_in = cellpose_in
            .combine(manifest_registration, by: 0)
            .map{ row -> tuple(
                row[0], // plate
                row[1], // well
                row[2], // row
                row[3], // col
                row[4], // nucl_channel
                row[5]  // cell_channel
            )}                
        } else {
            registration_out = null
        }
        
        // Run cellpose
        cellpose_out = cellpose(cellpose_in)
    
        //------------------------------------------------------------
        // Cellprofiler / get_features
        // TODO: BUG! Need to make sure the merge plate is done as well BEFORE adding into channel
        if (params.cpr_run) {
                        
            // Start with cellpose output
            // re-key channels
            cellpose_out = cellpose_out.map{row -> tuple(
                    row[0] + ":" + row[1], // key
                    row[0], // plate
                    row[1], // well
                    row[2], // row
                    row[3], // col
                    row[4], // cell masks
                    row[5]  // nucl masks
            )}
                
            //--------------------------------------------------------------------
            // Add registration
            if (registration_out != null) {
                // re-key output
                registration_out = registration_out.map{row -> tuple(
                    row[0] + ":" + row[1], // key
                    row[4], // merge plates
                    row[5], // path
                )} 
                
                // merge
                cellprofiler_in = cellpose_out.join(registration_out, by: 0)
                
            } else {
                cellprofiler_in = cellpose_out.map{row -> tuple(
                    row[0], // key
                    row[1], // plate
                    row[2], // well
                    row[3], // row
                    row[4], // col
                    row[5], // cell masks,
                    row[6], // nucl masks,
                    null,   // merge plates
                    file('NO_REGISTRATION'),   // registration path
                )}
                
            }
            
            //--------------------------------------------------------------------
            // Basicpy models        
            if (run_basicpy) {
            
                // Concat to plate string format
                basicpy_tmp = basicpy_out.map{ row -> tuple(
                    row[0], // plate
                    row[0] + "_ch" + row[1] + "=" + row[2] // plate : channel = path
                )}
    
                // debug only to test multiple channels
                //basicpy_out=basicpy_out.concat(basicpy_out)
    
                // Merge channel of same plates into single string 
                basicpy_tmp = basicpy_tmp
                .groupTuple(by:0)
                .map{row -> tuple(
                    row[1].join(" "), //  plate_ch1=path1 plate_ch2=path2 plate_chX=pathX
                    row[0] //plate
                )}
                
                // If there is a registration manifest, the basicpy output needs to be combined
                // into a single channel with the reference plate as the key
                if (params.rn_manifest_registration) {
                
                    // Read the manifest and turn it into a (qry, ref) channel
                    qry_ref=[]
                    new File(params.rn_manifest_registration).splitEachLine("\t") {fields ->
                        ref = fields[0]  
                        
                        // Add the ref plate
                        qry_ref.add([ref, ref])
                        
                        // Add the qry plates      
                        qry_items = fields[2].split(",")
                        for (qry in qry_items){
                            qry_ref.add([qry, ref])
                        }
                    }
                    
                    qry_ref_channel = Channel.from(qry_ref).map{row -> tuple(
                        row[1], // ref
                        row[0] // qry
                    )}
                    
                    // Add the reference plate to the basicpy output channel: (plate, bp_string, ref_plate)
                    basicpy_tmp = basicpy_tmp.combine(qry_ref_channel, by:1)
                    
                    // Group together by reference plate and concat the bp strings
                    basicpy_tmp = basicpy_tmp
                    .groupTuple(by:2)
                    .map{row -> tuple(
                        row[1].join(" "), //  plate_ch1=path1 plate_ch2=path2 plate_chX=pathX
                        row[2] //plate
                    )}
                }

                // Append the basicpy models for a plate into that channel              
                cellprofiler_in = cellprofiler_in.combine(basicpy_tmp, by: 1)
            } else {
                cellprofiler_in = cellprofiler_in.map{row -> tuple(
                        row[1], // plate
                        row[0], // key
                        row[2], // well
                        row[3], // row
                        row[4], // col
                        row[5], // cell masks
                        row[6], // nucl masks
                        row[7], // merge plates
                        row[8], // registration path
                        null    // basicpy models       
                )}
                
            }
            
            //--------------------------------------------------------------------
            // for in hybrid mode
            if (params.rn_hybrid) {
                cellprofiler_in = cellprofiler_in.combine(manifest, by: 0).map{row -> tuple(
                    row[0], // plate
                    row[1], // key
                    row[2], // well
                    row[3], // row
                    row[4], // col
                    row[5], // cell masks
                    row[6], // nucl masks
                    row[7], // merge plates
                    row[8], // registration path
                    row[9], // basicpy models   
                    (row[17] == "none") ? "none" : row[17].split(",").collect{it -> (it.toInteger() -1).toString()}.join(" ") // mask channels   
                )}
            } else {
                cellprofiler_in = cellprofiler_in.map{row -> tuple(
                    row[0], // plate
                    row[1], // key
                    row[2], // well
                    row[3], // row
                    row[4], // col
                    row[5], // cell masks
                    row[6], // nucl masks
                    row[7], // merge plates
                    row[8], // registration path
                    row[9], // basicpy models   
                    null // mask channels
                )}
            }
                
            cellprofiler_out = cellprofiler(cellprofiler_in, scaling_channel)
        
        }
        
}

