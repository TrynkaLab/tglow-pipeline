#!/bin/bash
 
set -e

# Load nextflow module
module load HGI/common/nextflow/23.10.0
 
#-------------------------------------------------------------------------------------
# General Nextflow variables #

export HTTP_PROXY='http://wwwcache.sanger.ac.uk:3128'
export HTTPS_PROXY='http://wwwcache.sanger.ac.uk:3128'
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms14G -Xmx14G -Dnxf.pool.maxThreads=2000"
 
# Nextflow version, change depending on installation you use
export NXF_VER=23.10.0
 
# Singularity cache folder to use (not currently used)
export NXF_SINGULARITY_CACHEDIR=/lustre/scratch123/hgi/mdt1/projects/healthy_imm_expr/resources/nextflow/cache/singularity
 
#-------------------------------------------------------------------------------------
# Path to nextflow file
NF_FILE=/software/teamtrynka/installs/tglow-pipeline/tglow.nf
 
# Arguments
# Output root folder
OUTPUT=""

# Plate level manifest
MANIFEST="manifest.tsv"

# Either lsf or lsf_inore_errors
PROFILE="lsf"

# Worflow name
WORKFLOW="$1"

# Run prefix for log
RUN_PREFIX=""

# Check if workflow was specified 
if [  -z "$OUTPUT" ]; then
    echo "[ERROR] No output root specified"
    exit 1
fi

# Check if workflow was specified 
if [  -z "$WORKFLOW" ]; then
    echo "[ERROR] No workflow specified"
    exit 1
fi

# Check if manifest exists
if [ ! -f "$MANIFEST" ]; then
    echo "[ERROR] Manifest file: ${MANIFEST} not found."
    exit 1
fi

# Check for blacklist, if not, create a dummy
if [ ! -f blacklist.tsv ]; then
    touch blacklist.tsv
fi

# Create output
mkdir -p ${OUTPUT}/workdir
mkdir -p ${OUTPUT}/results
mkdir -p ./logs

#-------------------------------------------------------------------------------------
# Build the command

# Set nextflow related arugments
CMD="nextflow \
-log logs/${RUN_PREFIX}${WORKFLOW}.nextflow.log  \
run ${NF_FILE} \
-profile ${PROFILE} \
-w ${OUTPUT}/workdir \
-resume \
-entry ${WORKFLOW}"

# Instance specific arguments go behind here
CMD="${CMD}  \
--rn_manifest ${MANIFEST} \
--rn_blacklist blacklist.tsv \
--rn_publish_dir ${OUTPUT}/results \
--cpr_pipeline_2d none \
--cpr_pipeline_3d none"

echo ${CMD}
eval ${CMD}