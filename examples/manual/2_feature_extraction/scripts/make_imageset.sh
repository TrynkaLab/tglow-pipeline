#!/bin/bash
# Usage
# Setup the batch file, change the arguments them
# ./make_imageset.sh

# Arguments
PROJECT=""
BATCH_NAME="template_pipeline"
BATCH="batch_files/${BATCH_NAME}.files"

# Main output folder to create subfolders in
OUT_DIR="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/${PROJECT}/2_feature_extraction/data/${BATCH_NAME}"

# Path with raw images
PROJECT_DIR="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/${PROJECT}"
IMG_RAW="${PROJECT_DIR}/1_pre_process/output/raw"

# Setup variables
PIPELINE_DIR="/software/teamtrynka/installs/tglow-core"
MERGER="${PIPELINE_DIR}/core/merge_tiffs_align.py"

#--------------------------------------------------------------------------------------
# Configure bash to work with conda
#source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
source /software/hgi/installs/conda-audited/miniconda/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/installs/tglow

while read file;
do
    WELL="$(echo $file | awk '{print $1}')"
    PLATE="$(echo $file | awk '{print $2}')"

    #echo $WELL $FOLDER
    
    mkdir -p ${OUT_DIR}/${PLATE}/${WELL}

    CMD="python ${MERGER} \
--well ${WELL} \
--input ${IMG_RAW} \
--output ${OUT_DIR}/ \
--plate ${PLATE} \
--max_project \
--no_zstack"

    echo "[INFO] Merging images:"
    echo "[INFO] ${CMD}"
    eval ${CMD}
    
done < $BATCH




