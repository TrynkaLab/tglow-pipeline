#!/usr/bin/env python 

import os
import re
import tifffile
import argparse
import numpy as np

from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
import logging

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Main loop
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Max project a plate/row/col/fields.ome.tiff into a YX tiff compatbible with cellprofiler.")

    parser.add_argument('-w', '--well', help='Well ID to merge', required=True)
    parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    parser.add_argument('-o','--output', help='Output prefix', required=True)
    parser.add_argument('-p','--plates', help='Subfolder in raw dir to process', nargs='+', required=True)
    parser.add_argument('-c', '--channel', help="The channel to save, defaults to the first", default=0)
    parser.add_argument('--suffix', help='Text to add at the end of the file, including extension', default=".tiff")
    parser.add_argument('--pattern', help="The pattern to index fields.", default="*.ome.tiff")
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    args = parser.parse_args()
        
    log.info(f"Input plates:\t{str(args.plates)}")
    log.info(f"Input fields:\t{str(args.fields)}")
    log.info(f"Input channel:\t{str(args.well)}")     
    log.info(f"Input well:\t{str(args.channel)}")     
    log.info(f"Input:\t\t{args.input}")
    log.info(f"Output:\t\t{args.output}")    
    log.info(f"Suffix:\t\t{args.suffix}")    
    log.info(f"Pattern:\t{args.pattern}")    

    reader = AICSImageReader(args.input, args.plates, fields_filter=args.fields, pattern=args.pattern)
    
    for plate in args.plates:
        for field in reader.fields.get(plate):
            log.info(f"Running field {field}")
            row, col = ImageQuery.well_id_to_index(args.well)
            iq = ImageQuery(plate, row, col, field, channel=int(args.channel))
            img = reader.read_image(iq)
            img = np.max(img, axis=0)
            
            outdir = f"{args.output}/{plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}"
            
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                log.info(f"Folder created: {outdir}")
                
            cur_out = f"{outdir}/{field}{args.suffix}"      
            tifffile.imwrite(cur_out, img, shape=img.shape, metadata={'axes': 'YX'}, compression="zlib", photometric='MINISBLACK')        

