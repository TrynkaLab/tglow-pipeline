#!/usr/bin/env python 

import numpy as np
import time
import logging
import argparse
import pandas as pd
import tifffile
import os
from tglow.io.tglow_io import AICSImageReader, AICSImageWriter
from tglow.io.image_query import ImageQuery

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Simple utility save a compressed version of a folder of existing images")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('-o','--output', help='Output folder', default="./")
    parser.add_argument('-p','--plates', help='Plates to process (at least one)', nargs='+', default=None)
    parser.add_argument('--pattern', help="The pattern to discover file,<pattern>.tiff", default="*_nucl_mask_*_cp_masks.tiff")
    args = parser.parse_args()

    #args.input = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/masks"
    #args.output = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/masks_compressed"
    #args.plates = "mo13_240518_TGlow_drugperturb_72h_plate1_cycle1"
    #args.pattern = "*_nucl_mask_*_cp_masks.tiff"
    
    reader = AICSImageReader(args.input, plates_filter=args.plates, pattern=args.pattern)
    
    # Some code to sanity check the new and old matrix are the same    
    #iq = ImageQuery(args.plates, 5, 5, 9)
    #cur_stack = reader.read_stack(iq)
    #reader2 = AICSImageReader(args.output, plates_filter=args.plates, pattern=args.pattern)
    #tmp = reader2.read_stack(iq)
    #np.array_equal(tmp, cur_stack)
    
    log.warning("This script will only save the first channel in a ZYX stack")
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
        
    log.info(f"Running {args.plates} with pattern {args.pattern}")
    
    for plate in reader.index.keys():
        cur_plate = reader.index[plate]
        for row in cur_plate.keys():
            cur_row = cur_plate[row]
            for col in cur_row.keys():
                cur_col = cur_row[col]
                
                log.info(f"Processing  plate {plate} row {ImageQuery.ID_TO_ROW[row]} col {col} ")
                for iq in cur_col.values():
                    cur_stack = reader.read_stack(iq)
                    cur_stack=cur_stack[0,]
                    
                    outdir = f"{args.output}/{iq.plate}/{iq.get_row_letter()}/{iq.col}/"
                    if not os.path.exists(outdir):
                        os.makedirs(outdir)
                        log.info(f"Folder created: {outdir}")
                                    
                    #writer.write_stack(cur_stack, iq)
                    tifffile.imwrite(f"{outdir}/{iq.field}{reader.suffix}", cur_stack, compression="zlib",  metadata={'axes': 'ZYX'})
    
    log.info("Successfully completed")