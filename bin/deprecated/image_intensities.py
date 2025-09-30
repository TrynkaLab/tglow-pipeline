#!/usr/bin/env python 

import numpy as np
import time
import logging
import argparse
import pandas as pd
import glob
import os
import warnings
import string

from skimage.measure import regionprops_table
from skimage.filters import threshold_otsu
from tglow.io.processed_image_provider import ProcessedImageProvider
from tglow.io.tglow_io import BlacklistReader, AICSImageReader, ControllistReader, ControlRecord
from tglow.io.image_query import ImageQuery
from tqdm import tqdm

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#---------------------------------------------------
# Extra properties for regionprops
def intensity_mean_nan(mask, intensity):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(intensity[mask])

def q1(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 1)

def q10(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 10)

def q25(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 25)

def q50(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 50)

def q75(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 75)

def q90(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 90)

def q99(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 99)

def q999(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 99.9)

def q9999(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 99.99)

def q99999(regionmask, intensity_image):
    return np.percentile(intensity_image[regionmask], 99.999)
#---------------------------------------------------


class ObjectIntensityCalculator():
    
    def __init__(self, input, obj_mask_dir, obj_mask_pattern, output, max_project, plates=None, plate_merge=None, flatfields=None, registration_dir=None, mask_channels=None, mask_dir=None,mask_pattern=None, blacklist=None, uint32=False):
        
        self.provider = ProcessedImageProvider(path=input,
                                               blacklist=blacklist,
                                               plate=plates,
                                               plate_merge=plate_merge,
                                               registration_dir=registration_dir,
                                               flatfields=flatfields,
                                               scaling_factors=None,
                                               mask_channels=mask_channels,
                                               mask_dir=mask_dir,
                                               mask_pattern=mask_pattern,
                                               uint32=uint32)
                    
        self.cell_mask_reader = AICSImageReader(path=obj_mask_dir, plates_filter=plates, blacklist=blacklist, pattern=obj_mask_pattern)
        self.max_project = max_project
        self.output=output
   
    def fetch_object_intensity(self, plate,  well):
        
        intensity_summary = None
        row, col = ImageQuery.well_id_to_index(well)
        fields = self.provider.plate_reader.fields_global[plate]
        pb = tqdm(total=len(fields), desc='Reading control images', unit='image')
        
        
        df = pd.DataFrame(colnames=["plate", "row", "col", "image", "channel", "q0", "q0.1", "q1", "q5", "q25", "q50", "q75", "q95", "q99", "q99.9", "q99.99", "q99.999", "q99.9999", "99.99999", "q100", "mean"])
        i = 0

        for field in self.provider.plate_reader.fields[plate][str(row)][str(col)]:
            
            iq = ImageQuery(plate, row, col, field)
            img = self.provider.fetch_image(iq)
            mask = self.cell_mask_reader.read_image(iq)
            pb.update(1)

            # Max project
            if self.max_project:
                img = np.max(img, axis=1)
                mask = np.max(mask, axis=1)[0,:,]
            else:
                mask = mask[0,]
            
            #extra_props = [intensity_mean_nan, q1, q10, q25, q50, q75, q90, q99, q999, q9999, q99999]    
            
            for channel in range(img.shape[0]):
                
                cur_img = img[channel,]
                cur_img[mask < 0 ] = np.nan
                df.loc[i,["q0", "q0.1", "q1", "q5", "q25", "q50", "q75", "q95", "q99", "q99.9", "q99.99", "q99.999", "q99.9999", "99.99999", "q100"]] = np.nanpercentile(cur_img,[0, 0.1, 1, 5, 25, 5, 75, 95, 99, 99.9, 99.99, 99.999, 99.9999, 99.99999, 100]).tolist()
                
                # Thresholded object intensity
                cur_img = img[channel,].astype(np.float32)
                thresh = threshold_otsu(cur_img)                

                
                #log.info(f"Identified {df.shape[0]} objects")
                df.loc[i,'plate'] = plate
                df.loc[i,'well'] = well
                df.loc[i,'field'] = field
                df.loc[i,'channel'] = channel

                if intensity_summary is None:
                    intensity_summary = df
                else:
                    intensity_summary = pd.concat([intensity_summary, df], axis=0)

        pb.close()
        
        return(intensity_summary)
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate per plate & per channel scaling factors based on intensity in foreground region")
    
    # Common options
    parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    parser.add_argument('-o','--output', help='Output prefix', required=True)
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+', required=True)
    parser.add_argument('-w','--well', help='Well to process', required=True)
    
    parser.add_argument('-m','--plate_merge', help='Plates to combine as multiple cycles of the same plate. Output will be seqeuntially added as additional channels in the order provided', nargs='+', default=None)
    parser.add_argument('--registration_dir', help="Path to registration root storing <plate>/<row>/<col>/<field>.pickle", default=None)
    parser.add_argument('--mask_dir', help="Path to mask root storing masks <plate>/<row>/<col>/<field><mask_pattern>.tiff", default=None)
    parser.add_argument('--mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py nucleus. <pattern>.tiff", default="*_nucl_mask_*_cp_masks.tiff")
    parser.add_argument('--mask_channels', help='Channels to mask, output will get 2 extra channels per mask_channel, in the order provided. This is applied after registering and flatfield correction, so use the final channel ids', nargs='+', default=None)
    parser.add_argument('--flatfields', help='Basicpy model dir for a channel. If merging channels ids are assigned seqeuntially for extra channels <plate>_ch<channel>=/path/to/model ', nargs='+', default=None)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint after applying basicpy model", action='store_true', default=False)

    # Specific options
    parser.add_argument('--obj_mask_dir', help='Base dir to cell masks organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('--obj_mask_pattern', help="The pattern to object discover masks, defaults to match run_cellpose.py cell. <pattern>.tiff", default="*_cell_mask_*_cp_masks.tiff")
    parser.add_argument('--max_project', help="Calculate flatfields on max projections. Automatically activates --all_planes", action='store_true', default=False)

    args = parser.parse_args()
    
    for plate in args.plate:
        calculator = ObjectIntensityCalculator(
                        input=args.input,
                        plates=[plate],
                        plate_merge=args.plate_merge,
                        output=f"{args.output}/{plate}",
                        blacklist=args.blacklist,
                        mask_dir=args.mask_dir,
                        mask_pattern=args.mask_pattern,
                        mask_channels=args.mask_channels,
                        obj_mask_dir=args.obj_mask_dir,
                        obj_mask_pattern=args.obj_mask_pattern,
                        registration_dir=args.registration_dir,
                        flatfields=args.flatfields,
                        uint32=args.uint32,
                        max_project=args.max_project)
        
        row, col = ImageQuery.well_id_to_index(args.well)
        row = f"{string.ascii_uppercase[int(row)-1]}"    
        
        result = calculator.fetch_object_intensity(plate, args.well)
        
        if not os.path.exists(f"{args.output}/{plate}/{row}/{col}/"):
            os.makedirs(f"{args.output}/{plate}/{plate}/{row}/{col}/")
            log.info(f"Folder created: {args.output}/{plate}/{row}/{col}/")
        
        
        result.to_csv(f"{args.output}/{plate}/{row}/{col}/{args.well}_object_intensities.tsv")

    
        result_filtered = result[result['threshold'] == 0]
        result_filtered.groupby('field')