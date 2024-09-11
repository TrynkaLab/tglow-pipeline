#!/bin/bash

# If something fails, exit
set -e

# LSF parameters
JS_GROUP="teamtrynka"
JS_MEM=16000
JS_CORES=1
JS_QUEUE="normal"
JS_TIME="12:00"

# Make a dir for the logs if it does not exist
mkdir -p logs

OUTPUT=../output/basicpy_models_all_plates

PLATES=(\
"<plate_name_1>" \
"<plate_name_2>")

# If you want to model per plate, uncomment this and fill the plate array above
plate=""
#for plate in ${PLATES[@]};
#do

for channel in 1 2 3 4 5;
do
    # Bsub command
    CMD="bsub -n ${JS_CORES} \
    -G ${JS_GROUP} \
    -M ${JS_MEM} \
    -q ${JS_QUEUE} \
    -W ${JS_TIME} \
    -R 'span[hosts=1] select[mem>${JS_MEM}] rusage[mem=${JS_MEM}]' \
    -o logs/ch${channel}-%J.out \
    -e logs/ch${channel}-%J.err \
    ./train_basicpy.sh ${OUTPUT}/${plate} ${channel} ${plate}"
    
    #echo $CMD
    # Execute the command
    eval ${CMD}
done
#done
