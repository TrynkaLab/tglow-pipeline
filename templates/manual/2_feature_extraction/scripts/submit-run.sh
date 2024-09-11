#!/bin/bash
# Usage
# submit-run.sh <run script.sh> <number of jobs in batch file>

# If something fails, exit
set -e

# LSF parameters
JS_GROUP="teamtrynka"
JS_MEM=6000
JS_CORES=1
JS_QUEUE="normal"
JS_TIME="2:00"
JS_NAME="$1"

# Make a dir for the logs if it does not exist
mkdir -p logs/${JS_NAME}

# Bsub command
CMD="bsub -n ${JS_CORES} \
-G ${JS_GROUP} \
-M ${JS_MEM} \
-q ${JS_QUEUE} \
-W ${JS_TIME} \
-R 'span[hosts=1] select[mem>${JS_MEM}] rusage[mem=${JS_MEM}]' \
-J '${JS_NAME}[1-${2}]%15' \
-o logs/${JS_NAME}/${JS_NAME}-%J-%I.out \
-e logs/${JS_NAME}/${JS_NAME}-%J-%I.err \
./${1}"

# Execute the command
#eval ${CMD}
echo ${CMD}