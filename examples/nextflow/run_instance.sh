#!/bin/bash
 
set -e
 
module load HGI/common/nextflow/23.10.0
 
export HTTP_PROXY='http://wwwcache.sanger.ac.uk:3128'
export HTTPS_PROXY='http://wwwcache.sanger.ac.uk:3128'
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms14G -Xmx14G -Dnxf.pool.maxThreads=2000"
 
# General Nextflow variables #
 
# Nextflow version, change depending on installation you use
export NXF_VER=23.10.0
 
# Singularity cache folder to use (not currently used)
export NXF_SINGULARITY_CACHEDIR=/lustre/scratch123/hgi/mdt1/projects/healthy_imm_expr/resources/nextflow/cache/singularity
 
# Path to nextflow file
NF_FILE=/software/teamtrynka/installs/tglow-core/nextflow/tglow.nf
 
# Output
OUTPUT="<OUTPUT>"
 
MANIFEST="<MANIFEST>"

WORKFLOW="$1"

mkdir -p ${OUTPUT}/workdir
mkdir -p ${OUTPUT}/results
 
# Build the command
# This first part should be the same for pretty much most pipelines
CMD="nextflow run ${NF_FILE} \
-profile lsf \
-w ${OUTPUT}/workdir \
-resume \
-entry ${WORKFLOW}"

# Pipeline specific arguments go behind here
CMD="${CMD}  \
--rn_manifest ${MANIFEST} \
--rn_publish_dir ${OUTPUT}/results \
--rn_image_dir ${OUTPUT}/results/images"

echo ${CMD}
eval ${CMD}