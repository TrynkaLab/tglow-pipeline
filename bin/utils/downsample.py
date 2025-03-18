#!/usr/bin/env python 

import os
import re
import tifffile
import argparse
import numpy as np

from aicsimageio.types import PhysicalPixelSizes
from skimage.transform import downscale_local_mean, resize
from tglow.io.tglow_io import AICSImageReader, AICSImageWriter
from tglow.io.image_query import ImageQuery
import logging

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Main loop
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Downsample a plate/row/col/fields.ome.tiff into a YX tiff ")

    parser.add_argument('-w', '--well', help='Well ID to merge', required=False, default=None)
    parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    parser.add_argument('-o','--output', help='Output prefix', required=True)
    parser.add_argument('-p','--plates', help='Subfolder in raw dir to process', nargs='+', required=True)
    parser.add_argument('-c', '--channel', help="The channel to save, defaults to the first", default=0)
    parser.add_argument('-z', '--downsample_z', help="The factor to downsample in z", default=1)
    parser.add_argument('-xy', '--downsample_xy', help="The factor to downsample in xy", default=1)
    parser.add_argument('--max_project', help="Output max projection over Z", action='store_true', default=False)
    parser.add_argument('--suffix', help='Text to add at the end of the file, including extension', default=".tiff")
    parser.add_argument('--pattern', help="The pattern to index fields.", default="*.ome.tiff")
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    parser.add_argument('--start_at_z', help="The first plane to use, zero indexed", default=None)
    parser.add_argument('--end_at_z', help="The last plane to use, zero indexed", default=None)
    args = parser.parse_args()
        
    log.info(f"Input plates:\t{str(args.plates)}")
    log.info(f"Input fields:\t{str(args.fields)}")
    log.info(f"Input channel:\t{str(args.channel)}")     
    log.info(f"Input well:\t\t{str(args.well)}")     
    log.info(f"Input:\t\t{args.input}")
    log.info(f"Output:\t\t{args.output}")    
    log.info(f"Suffix:\t\t{args.suffix}")    
    log.info(f"Pattern:\t\t{args.pattern}")    
    log.info(f"Downsample z:\t{args.downsample_z}")    
    log.info(f"Downsample xy:\t{args.downsample_xy}")    
    log.info(f"Max project:\t\t{args.max_project}")    
    log.info(f"Start plane:\t\t{args.start_at_z}")    
    log.info(f"End plane:\t\t{args.end_at_z}")   
    #args.input = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/Z_SPACING_OPT/pipeline/results/images"
    #args.plates = ["jm52_20250127_OPTMISATION_SPACING"]
    #args.fields = None
    #args.output =  "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/Z_SPACING_OPT/pipeline/results/images_7x"
    #args.pattern = "*.ome.tiff"
    #args.downsample_z = 7
    #args.downsample_xy = 0

    reader = AICSImageReader(args.input, args.plates, fields_filter=args.fields, pattern=args.pattern)
    writer = AICSImageWriter(args.output)

    downsample_z = int(args.downsample_z)
    downsample_xy = int(args.downsample_xy)

    if (args.start_at_z is not None):
        start_at_z = int(args.start_at_z)
    else:
        start_at_z = None

    if (args.end_at_z is not None):
        end_at_z = int(args.end_at_z)
    else:
        end_at_z = None


    for plate in args.plates:
        wells = []
        if args.well is not None:
            wells.append(args.well)
        else:
            wells = [item for item in reader.wells[plate]]
            
        for well in wells:
            row, col = ImageQuery.well_id_to_index(well)
            for field in reader.fields.get(plate):
                log.info(f"Running well {well} field {field}")
                
                iq = ImageQuery(plate, row, col, field)
                img = reader.read_image(iq)
                meta = reader.get_img(iq)
    
                log.info(f"Read image of shape: {img.shape} with dims {meta.dims}")
                
                if (start_at_z is not None and end_at_z is None):
                    img = img[:,range(start_at_z, img.shape[1]),:,:]
                    log.info(f"Dropped planes in range {start_at_z} to {end_at_z} to stack: {img.shape} with dims {meta.dims}")
                elif (start_at_z is None and end_at_z is not None):
                    img = img[:,range(0, end_at_z),:,:]
                    log.info(f"Dropped planes in range {start_at_z} to {end_at_z} to stack: {img.shape} with dims {meta.dims}")
                elif (start_at_z is not None and end_at_z is not None):
                    img = img[:,range(start_at_z, end_at_z),:,:]
                    log.info(f"Dropped planes in range {start_at_z} to {end_at_z} to stack: {img.shape} with dims {meta.dims}")
                
                pps = meta.physical_pixel_sizes
                
                # Downsample in z-direction
                if (downsample_z > 1):
                    planes_to_keep=[0]
                    for i in range(1, round(img.shape[1] / downsample_z)):
                        planes_to_keep.append((i * downsample_z)-1)

                    img = img[:,planes_to_keep,:,:]
                    log.info(f"Downsampled image {well} Z to {img.shape}")
                    pps = PhysicalPixelSizes(Z=pps.Z * downsample_z, X=pps.X, Y=pps.Y)

                # Downsample in xy-direction
                if (downsample_xy > 1):
                    img = downscale_local_mean(img, (1, 1, downsample_xy, downsample_xy))
                    log.info(f"Downsampled image {well} XY to {img.shape}")
                    pps = PhysicalPixelSizes(Z=pps.Z, X=pps.X * downsample_xy, Y=pps.Y * downsample_xy)

                # Max project
                if (args.max_project):
                    log.info("Max projecting")
                    img = np.max(img, axis=1, keepdims=True)
                
                log.info("Writing image stack")
                writer.write_stack(img, iq, channel_names=meta.channel_names, physical_pixel_sizes=pps)
               # outdir = f"{args.output}/{plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}"
               
            writer.write_image_stats(ImageQuery(plate, row, col, 1))
                
                #if not os.path.exists(outdir):
                #    os.makedirs(outdir)
                #    log.info(f"Folder created: {outdir}")
                    
                #cur_out = f"{outdir}/{field}{args.suffix}"      
                #tifffile.imwrite(cur_out, img, shape=img.shape, imagej=True, metadata={'axes': 'YX'})       
            
             

