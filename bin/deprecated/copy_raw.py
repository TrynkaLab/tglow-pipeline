#!/usr/bin/env python 

# Copied Martin Prete's scripts, edited for use here by TT
# Also edited the job names to make them easier to parse
####################################################################
# Script invocation:
# python actions/copy_raw.py \
#       --input_file="/path/to/raw/export/Images/Index.xml" \
#       --output_path="/lustre/scratch/raw" 
####################################################################
#
# - The script will bundle all the raw tiles that belong to each well
#   inside a file of filenames (fofn) and then submit a jobs to rsync
#   the contents of that file to the destination on lustre.
# - Why not copy as is? lustre doesn't like folders with a lot of
#   files in them. Each Phenix export can havve between 10-50k files.
#   That's why we create per-well folders to separate the raw tiles. 
# - The resulting Index.xml also will be modified o the images point
#   to the subfolders.
####################################################################

import os
import re
import string
import argparse
from xml.etree import ElementTree as ET
import time

def main(params):
    
    print(f"[+] Reading {params.input_file}")
    # read the index XML, we neeed to parse the namespace to read it with
    # python' ElementTree more easily	
    images_folder = os.path.dirname(params.input_file)
    metadata = ET.parse(params.input_file)
    default_namespace = re.findall(r'^{(.*)}',metadata.getroot().tag)[0]
    NS = {"PE": default_namespace}
    
		# collect plate information
    print(f"[+] Reading plate")
    plates = metadata.findall("./PE:Plates/PE:Plate",NS)
    assert len(plates)==1, "Expected only one plate"
    plate = plates[0]
    plate_name = plate.find("./PE:Name",NS).text
    plate_rows = int(plate.find("./PE:PlateRows",NS).text)
    plate_columns = int(plate.find("./PE:PlateColumns",NS).text)

    # gather all images
    images = {}
    # iterate over the images and group them by well
    print("[+] Grouping images by well")
    for image in metadata.findall("./PE:Images/PE:Image",NS):
        row = int(image.find("./PE:Row",NS).text)
        # convert numerical row to letter
        row = string.ascii_uppercase[row-1]
        column = image.find("./PE:Col",NS).text.zfill(2)
        well_key = f"{row}{column}"
        image_name = image.find('./PE:URL',NS)
        if well_key in images:
          images[well_key].append(image_name.text)
        else:
          images[well_key] = [image_name.text]
        # update image path inside Index.xml to refeerence
        # the image inside it's new path in the well-named folder
        image_name.text = os.path.join(well_key,image_name.text)
    print(f"[+] Wells = {len(images)}")
    
		# create files of filenames to use with rsync --files-from
    os.makedirs(os.path.join(params.output_path,plate_name,"fofn"), exist_ok=True)
    print("[+] Writing fofn and submitting rsync job")
    for k in images.keys():
      fofn = os.path.join(params.output_path,plate_name,"fofn",f"{k}.fofn")
      logs = os.path.join(params.output_path,plate_name,"logs")
      os.makedirs(logs, exist_ok=True) 
      with open(fofn, "wt") as f:
        f.write("\n".join(images[k]))
      # submit a job per well to copy in parallel as much as possible
      os.system(f'bsub -G teamtrynka -J {plate_name}_well{k} -q imaging -n 1 -M 50M -R "select[mem>50M] rusage[mem=50M]" -o "{logs}/copy_raw-%J.log" -e "{logs}/copy_raw-%J.log" rsync -arvh --files-from="{fofn}" "{images_folder}" "{os.path.join(params.output_path,plate_name,k)}"')
      # Sleep for 4 seconds in between jobs not to overload the lustre
      time.sleep(2)
      
    print("[+] Writing update Index.xml")
    idx = os.path.join(params.output_path,plate_name,"Index.xml")
    ET.register_namespace('', default_namespace)
    ET.indent(metadata)
    metadata.write(idx, encoding="utf-8")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Copy raw tiles from Perkin Elmer export to individual well folders')
    parser.add_argument('--input_file', type=str, required=False, help='Path to the PE Index file')
    parser.add_argument('--output_path', type=str, required=True, help='Path to the output directory where to copy the files and creaate the folder structure')
    try:
       args = parser.parse_args()
    except:
       parser.print_help()
       exit(1)

    main(args)

