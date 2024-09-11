#!/bin/bash

# If something fails, exit
set -e

# LSF parameters
JS_GROUP="teamtrynka"
JS_MEM=4000
JS_CORES=1
JS_QUEUE="normal"

JS_TIME="12:00"

# Make a dir for the logs if it does not exist
mkdir -p logs

# Bsub command
CMD="bsub -n ${JS_CORES} \
-G ${JS_GROUP} \
-M ${JS_MEM} \
-q ${JS_QUEUE} \
-W ${JS_TIME} \
-R 'span[hosts=1] select[mem>${JS_MEM}] rusage[mem=${JS_MEM}]' \
-o logs/tglow-nextflow-main-%J.out \
-e logs/tglow-nextflow-main-%J.err \
./run_instance.sh"

# Execute the command
eval ${CMD}
