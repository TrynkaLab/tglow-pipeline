#!/bin/bash

# If something fails, exit
set -e

#------------------------------------------------------------------------------------
# Project settings
PROJECT="<PROJECT>"
PROJECT_DIR="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/${PROJECT}/"

#  Main output folder to create subfolders in
OUT_DIR="${PROJECT_DIR}/2_feature_extraction/output/<OUTPUT>"

# Job array file, one job per line <well> <plate>
ARRAY_FILE="${PROJECT_DIR}/2_feature_extraction/scripts/batch_files/<BATCH FILE>"

# Software variables 
PIPELINE_DIR="/software/teamtrynka/installs/tglow-core"
PIPELINE="${PIPELINE_DIR}/cellprofiler/${PROJECT}/<.CPPIPE>"

# Path with raw images
IMG_RAW="${PROJECT_DIR}/1_pre_process/output/raw"

#------------------------------------------------------------------------------------
# Software settings
MERGER="${PIPELINE_DIR}/core/merge_tiffs_align.py"
CP_PLUGINS="/software/teamtrynka/installs/cellprofiler/plugins"
set _JAVA_OPTIONS=-Xmx7g

#------------------------------------------------------------------------------------
# Well and plate (subfolder) to run
# Grab from jobindex
WELL="$(sed -n "${LSB_JOBINDEX}p" ${ARRAY_FILE} | awk '{print $1}')"
PLATE="$(sed -n "${LSB_JOBINDEX}p" ${ARRAY_FILE} | awk '{print $2}')"
PLATE_MERGE="$(sed -n "${LSB_JOBINDEX}p" ${ARRAY_FILE} | awk '{print $3}')"

#------------------------------------------------------------------------------------
# Arguments
PLANES="1 2 3 4 5 6 7 8 9 10 11 12 13 14"
FIELDS="1 2 3 4 5 6 7 8 9"
CHANNELS="1 2 3 4 5"

# Path with pretrained basicpy models
BASICPY_DIR="${PROJECT_DIR}/1_pre_process/output/basicpy_models"

# Supply basicpy models, if used comment blank below, if no models to be used, use blank below
BASICPY="--basicpy_model \
5=${BASICPY_DIR}/<MODEL_SUBDIR> \
8=${BASICPY_DIR}/<MODEL_SUBDIR>"

BASICPY=""

# If aligning multiple cycles comment blank below, otherwise use blank below
MERGE_ARGS="--plate_merge ${PLATE_MERGE} \
--ref_channel 5 \
--qry_channel 3 \
--channels_merge 1 2 3 4"

MERGE_ARGS=""

# Should merged images be removed after cellprofiler is done
CLEAN=true

# [CAREFULL!] Should the previous output be wiped BEFORE running
WIPE=false


#------------------------------------------------------------------------------------
# Start the run
echo "[INFO] -------------------------------------------------------------"
echo "[INFO] Well: ${WELL} Plate: ${PLATE}"

if [ -z "${WELL}" ]; then
    echo "[ERROR] No well id supplied"
    exit 1
fi

if [ -z "${PLATE}" ]; then
    echo "[ERROR] No subfolder / plate supplied"
    exit 1
fi

# Output dir
OUT_DIR=${OUT_DIR}/${PLATE}/${WELL}

# Wipe previous output
if [ "${WIPE}" = true ]; then
    if [ -d "${OUT_DIR}" ]; then
        echo "[INFO] Wiping previous run output"
        rm -r ${OUT_DIR}
    fi
fi

# Check if there is a previously completed run
if [ -f  "${OUT_DIR}/.done" ]; then
    echo "[INFO] .done file found, skipping job. If you want to run anyway, set WIPE=true"
    exit 0;
fi

# Configure bash to work with conda
#source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
source /software/hgi/installs/conda-audited/miniconda/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/installs/tglow

# Check if input exists
if [ -d "${IMG_RAW}/${PLATE}/${WELL}" ]; then
    echo "[INFO] Input folder check passed."
else
    echo "[ERROR] Raw image input folder not detected, exiting."
    echo "[ERROR] ${IMG_RAW}/${PLATE}/${WELL}"
    exit 1
fi

# Dir to store merged images
IMG_DIR=${OUT_DIR}/images

echo "[INFO] -------------------------------------------------------------"

if [ -d "${IMG_DIR}" ]; then
    echo "[INFO] Image folder detected. Skipping merging of images"
else

    # Do some stuff to fetch and merge images for well
    mkdir -p ${IMG_DIR}
    
    # Command for merging
    CMD_MERGE="python ${MERGER} \
    --well ${WELL} \
    --input ${IMG_RAW} \
    --output ${OUT_DIR}/images \
    --plate ${PLATE} \
    --planes ${PLANES} \
    --fields ${FIELDS} \
    --channels ${CHANNELS} \
    --no_zstack \
    --max_project \
    ${BASICPY} \
    ${MERGE_ARGS}"

    echo "[INFO] Merging images:"
    echo "[INFO] ${CMD_MERGE}"
    eval ${CMD_MERGE}
fi

echo "[INFO] -------------------------------------------------------------"

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/installs/cellprofiler

# Run cell profiler
CMD="cellprofiler \
-c \
-r \
-p ${PIPELINE} \
-o ${OUT_DIR} \
-i ${IMG_DIR} \
--plugins-directory ${CP_PLUGINS}"
    
echo "[INFO] Running CellProfiler:"
echo "[INFO] ${CMD}"
eval ${CMD}
echo "[INFO] -------------------------------------------------------------"

if [ "${CLEAN}" = true ]; then
    echo "[INFO] Cleanup, removing merged images and logs"
    rm -r ${IMG_DIR}
    #rm "logs/${1}/*-${LSB_JOBINDEX}.out"
    #rm "logs/${1}/*-${LSB_JOBINDEX}.err"
fi

# File to confirm successfull completion
touch  "${OUT_DIR}/.done"