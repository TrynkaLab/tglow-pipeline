#!/bin/bash
# Usage
# run_companion_ome.sh <batch_file> <stitched input dir>
# batch file should be a list of plate names, one per line

# Configure bash to work with conda
#source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
source /software/hgi/installs/conda-audited/miniconda/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/installs/cellprofiler

PIPELINE_DIR="/software/teamtrynka/installs/tglow-core"
BATCH="batch_files/<batch>.files"
INPUT_DIR="../output/stitched"

while read batch;
do

path="${INPUT_DIR}/${batch}/index.xml"

CMD="python ${PIPELINE_DIR}/core/companion.py \
--input_file ${path}"

echo $CMD
eval $CMD

done < $BATCH
