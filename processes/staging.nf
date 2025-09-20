#!/usr/bin/env nextflow

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
    stub:
        """
        mkdir -p "$plate/$row/$col"
        cd "$plate/$row/$col"
        touch 1.ome.tiff
        touch CHECKSUMS.txt
        """         
}


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
        cmd
    
    stub:
        """
        touch manifest.tsv
        touch Index.json
        touch Index.xml
        touch acquisition_info.txt
        """    
}

// Create a manifest for a dir of images if it does not exist yet
process index_imagedir {
    label "normal"
    conda params.tg_conda_env
    //storeDir "${params.rn_publish_dir}/$input_dir/", saveAs: { filename -> filename.split('/')[-1] }
    //storeDir "${params.rn_publish_dir}/$input_dir"
    publishDir "${params.rn_publish_dir}/$input_dir", mode: "copy"
    input:
        val previous_completed
        val input_dir
        path images, stageAs: "input_images"
        val plate
    output:
        path "$plate/manifest.tsv"
    script:
        """
        index_folder.py \
        --input input_images \
        --plates $plate \
        --output ./
        """
    stub:
        """
        mkdir -p $plate
        cd $plate
        touch manifest.tsv
        """
}