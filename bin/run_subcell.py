#!/usr/bin/env python 

import os
import tifffile
import argparse
import numpy as np

import logging
from tglow.io.image_query import ImageQuery
from tglow.io.processed_image_provider import ProcessedImageProvider
from tglow.io.tglow_io import AICSImageWriter, AICSImageReader
from skimage.transform import rescale
# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

from skimage import measure
from skimage import util
import argparse
import datetime
import logging
import os
import sys
import pandas as pd
import requests
import torch
import yaml
from skimage.io import imread
import boto3
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError
from urllib.parse import urlparse

import inference
from vit_model import ViTPoolClassifier


os.environ["DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"


def get_crops(img, msk, outdim=640):
    regions = measure.regionprops(msk)
    
    crops = {}

    padding=outdim/2

    for i, region in enumerate(regions):
        z, y, x = region.centroid
    
        id = f"{int(y)}:{int(x)}"
    
        crop = np.zeros((img.shape[0], 640, 640))
    
        minr = int(max(0, y - padding))
        minc = int(max(0, x - padding))
        
        if minr == 0:
            padding_y = int(y)
        else:
            padding_y = padding
            
        if minc == 0:
            padding_x = int(x)
        else:
            padding_x=padding
        
        maxr = int(min(img.shape[1], y + padding_y))
        maxc = int(min(img.shape[2], x + padding_x))
        
        if maxr == img.shape[1]:
            minr == y - (maxr-y)
        
        if maxc == img.shape[1]:
            minc == x - (maxr-x)
        
        curd = img[:,minr:maxr, minc:maxc]
        
        large_center = np.array(crop.shape) // 2
        small_center = np.array(curd.shape) // 2
        
        start = large_center - small_center
        end = start + curd.shape
        
        crop[start[0]:end[0], start[1]:end[1], start[2]:end[2]] = curd
                
        # Crop the image
        #crops.append(crop)
        crops[id] = crop
        
    return(crops)
    
# Main loop
if __name__ == "__main__":
            
    # CLI 
    # parser.add_argument('-w', '--well', help='Well ID to merge', required=True)
    
    # # Common argulents
    # parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    # parser.add_argument('-o','--output', help='Output prefix', required=True)
    # parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+', required=True)
    # parser.add_argument('-w','--well', help='Subfolder in raw dir to process', nargs='+', required=True)
    # #parser.add_argument('--mask_dir', help="Path to mask root storing masks <plate>/<row>/<col>/<field><mask_pattern>.tiff", default=None)
    # parser.add_argument('--mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py cell. <pattern>.tiff", default="*_cell_mask_*_cp_masks.tiff")
    # parser.add_argument('--scale_factor', help="The scale factor to get to 80nm per pixel", default=149/80)
    #parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. <plane #1> | [<plane #1> <plane #2> <plane #n>]', nargs='+', default=None)
    
    #cd /lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/test_subcell
    parser = argparse.ArgumentParser(description="Run subcell inference on a plate/row/col/field.ome.tiff and associated masks ")
    args = parser.parse_args()
    args.output="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/test_subcell"
    args.input="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/processed_images"
    args.well = "C10"
    args.plate = "mo13_240518_TGlow_drugperturb_72h_plate1_cycle1"
    args.mask_pattern="*_cell_mask_*_cp_masks.tiff"
    args.scale_factor=149/80
    
    
    img = tifffile.imread(f"{args.input}/{args.plate}/C/10/1.ome.tiff")
 
    img_reader = AICSImageReader(args.input, plates_filter=args.plate)
    msk_reader = AICSImageReader(args.input, plates_filter=args.plate, pattern=args.mask_pattern)

    row, col = ImageQuery.well_id_to_index(args.well)
    iq = ImageQuery(args.plate, row, col, 1)
    
    # Assume images are already MP
    img = img_reader.read_stack(iq)[:,0,:,:]
    msk = msk_reader.read_stack(iq)[:,0,:,:]
    
    img = rescale(img, args.scale_factor, channel_axis=0, order=0)   
    msk = rescale(msk, args.scale_factor, channel_axis=0, order=0)   
    
    crops = get_crops(img, msk)    
    
    # This is the general configuration variable. We are going to use the special key "log" in the dictionary to use the log in our code
    config = {"log": log}

    # If you want to use constants with your script, add them here
    config["model_channels"] = "rybg"
    config["model_type"] = "mae_contrast_supcon_model"
    config["update_model"] = False
    config["create_csv"] = False
    config["gpu"] = -1

    with open(
        os.path.join(
            "models",
            config["model_channels"],
            config["model_type"],
            "model_config.yaml",
        ),
        "r",
    ) as config_buffer:
        model_config_file = yaml.safe_load(config_buffer)

    classifier_paths = None
    
    classifier_paths = None
    if "classifier_paths" in model_config_file:
        classifier_paths = model_config_file["classifier_paths"]
    encoder_path = model_config_file["encoder_path"]

    #--------------------------
    model_config = model_config_file.get("model_config")
    model = ViTPoolClassifier(model_config)
    model.load_model_dict(encoder_path, classifier_paths)
    model.eval()
      
    if torch.cuda.is_available() and config["gpu"] != -1:
        device = torch.device("cuda:" + str(config["gpu"]))
    else:
        config["log"].warning("CUDA not available. Using CPU.")
        device = torch.device("cpu")
    model.to(device)
  
  
    localization_channel = 1
  
  
    # Loop over images here
    cell_data=crops[4][(8,5,7,1), :,:]
    
    cell_data = [[cell_data[i]] for i in range(cell_data.shape[0])]
    #cell_crop = np.stack(cell_data, axis=1)

    
    embedding, probabilities = inference.run_model(
        model,
        cell_data,
        device,
        os.path.join("./testing_1"),
    )
    
    
    cell_crop = np.stack(cell_data, axis=1)
    cell_crop = torch.from_numpy(cell_crop).float().to(device)
    cell_crop = min_max_standardize(cell_crop)

    output = model(cell_crop)


    save_attention_map(output.attentions, (cell_crop.shape[2], cell_crop.shape[3]), "testing_1_12_attn")
    
        
    curr_probs_l = probabilities.tolist()
    max_location_class = curr_probs_l.index(max(curr_probs_l))
    max_location_name = inference.CLASS2NAME[max_location_class]
    max_3_location_classes = sorted(
        range(len(curr_probs_l)), key=lambda sub: curr_probs_l[sub]
    )[-3:]
    max_3_location_classes.reverse()
    max_3_location_names = (
        inference.CLASS2NAME[max_3_location_classes[0]]
        + ","
        + inference.CLASS2NAME[max_3_location_classes[1]]
        + ","
        + inference.CLASS2NAME[max_3_location_classes[2]]
    )
    
    
    
    
    
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    
    log.info("-----------------------------------------------------------")
    
    runner.run(args.well)

