#!/bin/sh
# Usage
# run_omero_upload.sh <RENDER_MANIFEST> <UPLOAD_MANIFEST>

# Configure bash to work with conda
software/hgi/installs/anaconda3/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/installs/omero/

RENDER_MANIFEST=$1
UPLOAD_MANIFEST=$2 

#RENDER_MANIFEST=batch_files/ki67_test_omero_manifest_full.csv
#UPLOAD_MANIFEST=/lustre/scratch125/humgen/projects/cell_activation_tc/projects/KI67_TEST/1_pre_process/scripts/batch_files/ki67_test_omero_manifest.txt
#omero_group="JM_TCA"
#omero_username="ob7"
#omero_project="ki67_test"

while read line;
do
file="$(echo $line | awk 'BEGIN{FS=","}{print $1}').render.yml"

echo $file

#file="test.yml"
echo -e "channels:" > $file

echo -e "  1:" >> $file
echo -e "    label: \"$(echo $line | awk 'BEGIN{FS=","}{print $7}')\"" >> $file
echo -e "    start: $(echo $line | awk 'BEGIN{FS=","}{print $8}')" >> $file
echo -e "    end: $(echo $line | awk 'BEGIN{FS=","}{print $9}')" >> $file
echo -e "    color: \"$(echo $line | awk 'BEGIN{FS=","}{print $10}')\"" >> $file

echo -e "  2:" >> $file
echo -e "    label: \"$(echo $line | awk 'BEGIN{FS=","}{print $11}')\"" >> $file
echo -e "    start: $(echo $line | awk 'BEGIN{FS=","}{print $12}')" >> $file
echo -e "    end: $(echo $line | awk 'BEGIN{FS=","}{print $13}')" >> $file
echo -e "    color: \"$(echo $line | awk 'BEGIN{FS=","}{print $14}')\"" >> $file

echo -e "  3:" >> $file
echo -e "    label: \"$(echo $line | awk 'BEGIN{FS=","}{print $15}')\"" >> $file
echo -e "    start: $(echo $line | awk 'BEGIN{FS=","}{print $16}')" >> $file
echo -e "    end: $(echo $line | awk 'BEGIN{FS=","}{print $17}')" >> $file
echo -e "    color: \"$(echo $line | awk 'BEGIN{FS=","}{print $18}')\"" >> $file

echo -e "  4:" >> $file
echo -e "    label: \"$(echo $line | awk 'BEGIN{FS=","}{print $19}')\"" >> $file
echo -e "    start: $(echo $line | awk 'BEGIN{FS=","}{print $20}')" >> $file
echo -e "    end: $(echo $line | awk 'BEGIN{FS=","}{print $21}')" >> $file
echo -e "    color: \"$(echo $line | awk 'BEGIN{FS=","}{print $22}')\"" >> $file

echo -e "  5:" >> $file
echo -e "    label: \"$(echo $line | awk 'BEGIN{FS=","}{print $23}')\"" >> $file
echo -e "    start: $(echo $line | awk 'BEGIN{FS=","}{print $24}')" >> $file
echo -e "    end: $(echo $line | awk 'BEGIN{FS=","}{print $25}')" >> $file
echo -e "    color: \"$(echo $line | awk 'BEGIN{FS=","}{print $26}')\"" >> $file

done < $RENDER_MANIFEST

cd /software/teamtrynka/installs/omero/src
python client.py --manifest $UPLOAD_MANIFEST

#python client.py \
#--file /lustre/scratch125/humgen/projects/cell_activation_tc/projects/KI67_TEST/1_pre_process/output/stitched/3012024_mo13_TGlow_ki67Test_30MinPerm/3012024_mo13_TGlow_ki67Test_30MinPerm.companion.ome \
#--omero_group JM_TCA --omero_project ki67_test --omero_dataset 3012024_mo13_TGlow_ki67Test_30MinPerm --omero_username ob7 

#[25/106] F2_F1.ome.tiff
#[+] validating
#[x] error: file 'F2_F1.ome.tiff' already exists in group 'JM_TCA'
#---------
#[26/106] F3_F1.ome.tiff
#[+] validating
#[x] error: file 'F3_F1.ome.tiff' already exists in group 'JM_TCA'
#---------
#[27/106] F4_F1.ome.tiff
#[+] validating
#[x] error: file 'F4_F1.ome.tiff' already exists in group 'JM_TCA'
#---------
#[28/106] F5_F1.ome.tiff
#[+] validating
#[x] error: file 'F5_F1.ome.tiff' already exists in group 'JM_TCA'
#---------
#[29/106] F6_F1.ome.tiff
#[+] validating
#[x] error: file 'F6_F1.ome.tiff' already exists in group 'JM_TCA'

#python client.py --file /lustre/scratch125/humgen/projects/cell_activation_tc/projects/KI67_TEST/1_pre_process/output/stitched/3012024_mo13_TGlow_ki67Test_10minperm/F4_F1.ome.tiff --omero_group JM_TCA --omero_project test --omero_dataset testset2 --omero_username ob7 