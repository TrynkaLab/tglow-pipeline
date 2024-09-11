#!/bin/bash
# Usage
# run_stitch.sh <MANIFEST>

set -e

module load cellgen/acapella

manifest=$1
manifest_entry="$(tail -n +2 ${manifest} | sed -n ${LSB_JOBINDEX}p)"

index_file="$(echo $manifest_entry | awk 'BEGIN {FS=","}{print $1}')"
out_folder="$(echo $manifest_entry | awk 'BEGIN {FS=","}{print $2}')"
z_projection="$(echo $manifest_entry | awk 'BEGIN {FS=","}{print $3}')"
gap="$(echo $manifest_entry | awk 'BEGIN {FS=","}{print $4}')"
planes="ALL"

mkdir -p "${out_folder}"

CMD="stitch.sh \
--index-file \"${index_file}\" \
--output-dir \"${out_folder}\" \
--z-projection ${z_projection} \
--silent false \
--wells ALL \
--gap ${gap} \
--planes ${planes} \
--channels ALL \
--fields ALL \
--timepoints ALL"

echo $CMD
eval $CMD