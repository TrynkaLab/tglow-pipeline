#!/usr/bin/env nextflow

process index_imagedir {
    label "tiny"
    conda params.tg_conda_env
    storeDir "${params.rn_publish_dir}/$input_dir/" 
    
    input:
        val input_dir
    output:
        path "manifests/*"
    script:
    """
    index_folder.py \
    --input ${params.rn_publish_dir}/$input_dir/ \
    --output ./manifests
    """
}


process fetch_model {
    label "midmem"
    conda params.sc_conda_env
    storeDir "${params.rn_publish_dir}/subcell/models" 
    input:
        val channels
        val type
    output:
        path "models/*"
    script:
        cmd = 
        """
        mkdir -p models/$channels/$type
        
        cp ${file(assets/subcell/models/$channels/$type/model_config.yaml)} ./models/$channels/$type/
        
        fetch_subcell_models.py \
        --input models/$channels/$type/model_config.yaml  \
        --model_urls ${file('assets/subcell/models_urls.yaml')}
        --model_channels $channels \
        --model_type $type \
        --output ./
        """
        cmd
}


process subcell {
    label "midmem"
    conda params.sc_conda_env
    storeDir "${params.rn_publish_dir}/subcell/results" 
    
    input:
        tuple val(well), val(row), val(col), val(plate), val(index_xml)
        path model_dir
        path input_dir
    output:
        path "$plate/$row/$col/*"
        
    script:
        cmd=
        """
        run_subcell.py \
        --input $input_dir \
        --model $model_dir \
        --well $well \
        --output ./features \
        
        """
        
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