#!/bin/bash
# Usage
# run_copy_raw.sh <batch_file> <output dir>
# batch file should be a list of plate names, one per line

# copy_raw.py authored by Martin Prete
# Script invocation:
# python actions/copy_raw.py \
#       --input_file="/path/to/raw/export/Images/Index.xml" \
#       --output_path="/lustre/scratch/raw" 

# Configure bash to work with conda
#source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
source /software/hgi/installs/conda-audited/miniconda/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/installs/tglow

# Setup variables
PIPELINE_DIR="/software/teamtrynka/installs/tglow-core"
BATCH="$1"
OUT_DIR="$2"

while read plate;
do
echo "[INFO] ${file}"

# For Operetta exports use Index.idx.xml here insteaf
python ${PIPELINE_DIR}/core/copy_raw.py \
--input_file="/nfs/t217_imaging/HarmonyExports/${plate}/Images/Index.xml" \
--output_path="${OUT_DIR}"

done < $BATCH
