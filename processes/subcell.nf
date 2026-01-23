#!/usr/bin/env nextflow


process fetch_model {
    label "normal"
    conda params.sc_dl_conda_env
    container params.sc_dl_container
    
    storeDir "${params.rn_publish_dir}/subcell/" 
    input:
        val channels
        val type
        path config
        path url
    output:
        path "models/$channels/$type"
    script:
        cmd = 
        """
        mkdir -p models/$channels/$type
        
        cp $config ./models/$channels/$type/
        
        fetch_subcell_models.py \
        --input models/$channels/$type/model_config.yaml  \
        --model_urls $url \
        --model_channels $channels \
        --model_type $type \
        --output ./
        """
        cmd
    stub:
        """
        mkdir -p models/$channels/$type
        cd models/$channels/$type
        touch model_config.yaml
        """
}


process subcell {
    label params.sc_label
    
    conda params.sc_conda_env
    container params.sc_container
    
    storeDir "${params.rn_publish_dir}/subcell/" 
    
    input:
        tuple val(well), val(row), val(col), val(plate), val(index_xml)
        path model_dir
        path input_dir
    output:
        path "features/$plate/$row/$col/*"
        
    script:
        cmd=
        """
        run_subcell.py \
        --input $input_dir \
        --model $model_dir \
        --plate $plate \
        --well $well \
        --output ./features \
        --channels $params.sc_channels \
        --ch_nucl $params.sc_nucl \
        --ch_tub $params.sc_tub \
        --ch_er $params.sc_er \
        --scale_factor $params.sc_scale \
        --mask_pattern '*_cell_mask_*_cp_masks.tiff'\
        """
        
        if (params.sc_gpu) {
            cmd += " --gpu"
        }
        
        if (params.sc_dont_mask) {
            cmd += " --dont_mask"
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
    stub:
        """
        mkdir -p features/$plate/$row/$col
        cd features/$plate/$row/$col
        touch ${plate}_${well}.txt
        """
}