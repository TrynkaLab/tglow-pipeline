#!/usr/bin/env python 

import numpy as np
import time
import logging
import argparse
import pandas as pd
import glob
import os
from tglow.io.tglow_io import BlacklistReader, AICSImageReader, ControllistReader, ControlRecord
from tqdm import tqdm

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ScalingCalculator():
    
    def __init__(self, controls, input, masks, output, max_project, mask_pattern,  plates=None, blacklist=None, fields=None):
        
        if blacklist is not None:
            blacklist = BlacklistReader(blacklist).read_blacklist()
        
        self.controls = ControllistReader(controls, plates_filter=plates, blacklist=blacklist).read_controlist()
        self.image_reader = AICSImageReader(path=input, plates_filter=plates, blacklist=blacklist)
        self.mask_reader = AICSImageReader(path=masks, plates_filter=plates, blacklist=blacklist, pattern=mask_pattern)
        self.max_project = max_project
        self.output=output

    def make_intensity_summary(self):
                
        if self.nobjects is None:
            self.fetch_object_intensity()
            
        nobjects=self.nobjects
        object_size=self.object_size
        intensity=self.intensity
        med_intensity=self.med_intensity
        max_intensity=self.max_intensity
        
        # Fetch unique plates
        plate_channel = set()
        for control in self.controls:
            for channel in control.channels:
                plate_channel.add(f"{control.plate}:{channel}")
        
        mean_intensity={}
        integrated_intensity={}
        
        df = pd.DataFrame(columns=["plate", "channel", "nobjects", 'nfields', "mean_mean_object_intensity", "mean_integrated_intensity", "mean_median_intensity", "mean_max_intensity"])
        
        i = 0
        for pc in plate_channel:
            plate = pc.split(":")[0]
            channel = int(pc.split(":")[1])
            
            # Average average intensity in volume per field
            mean_intensity[pc] = np.sum(np.array(intensity[plate][channel]) / np.array(object_size[plate][channel])) / len(intensity[plate][channel])
            
            # Average integrated intensity per field
            integrated_intensity[pc] = np.sum(np.array(intensity[plate][channel])) / len(intensity[plate][channel])
            mean_median_intensity = np.sum(np.array(med_intensity[plate][channel])) / len(med_intensity[plate][channel])
            mean_max_intensity = np.sum(np.array(max_intensity[plate][channel])) / len(max_intensity[plate][channel])

            
            df.at[i, "plate"] = plate
            df.at[i, "channel"] = channel
            df.at[i, "nobjects"] = np.sum(np.array(nobjects[plate][channel]))
            df.at[i, "nfields"] = len(np.array(nobjects[plate][channel]))
            df.at[i, "mean_mean_object_intensity"] = mean_intensity[pc]
            df.at[i, "mean_integrated_intensity"] = integrated_intensity[pc]
            df.at[i, "mean_median_intensity"] = mean_median_intensity
            df.at[i, "mean_max_intensity"] = mean_max_intensity

            log.info(f"Plate: {plate}, channel: {channel}")
            log.info(f"Mean mean intensity per field: {mean_intensity[pc]}")
            log.info(f"Mean integrated intensity per field: {integrated_intensity[pc]}")
            i+=1
            
        df.to_csv(f"{self.output}/plate_level_control_intensity_summary.tsv", sep='\t', index=False)

        #return mean_intensity, integrated_intensity
        

    def fetch_object_intensity(self):
        
        nobjects = {}
        object_size = {}
        intensity = {}
        max_intensity = {}
        med_intensity = {}

        pb = tqdm(total=len(self.controls * 9), desc='Reading control images', unit='image')
        for control in self.controls:
            
            for field in self.image_reader.fields[control.plate]:
                img = self.image_reader.read_image(control.get_query(field))
                mask = self.mask_reader.read_image(control.get_query(field))
                pb.update(1)

                # Keep track of the number of objects
                cur_objects=np.max(mask)

                # Binarize mask
                mask = mask > 0
                cur_area=np.sum(mask)

                # Max project
                if self.max_project:
                    img = np.max(img, axis=1)
                    mask = np.max(mask, axis=1)[0,:,]
                    #log.debug(f"Max projected img and mask to {img.shape} and {mask.shape}")
                else:
                    mask= mask[0,]
                
                for channel in control.channels:
                    # Set masked region in a channel to zero so non-object pixels are not
                    # considered for summing up the intensity                    
                    img[channel,][mask <=0] = 0 
                    
                    # Save the number of objects and the total object area
                    if control.plate not in nobjects.keys():
                        nobjects[control.plate] = {}
                        object_size[control.plate] = {}

                    if channel not in nobjects[control.plate].keys():
                            nobjects[control.plate][channel] = []
                            object_size[control.plate][channel] = []
                             
                    nobjects[control.plate][channel].append(cur_objects)
                    object_size[control.plate][channel].append(cur_area)

                    # Save the intensity
                    if control.plate not in intensity.keys():
                        intensity[control.plate] = {}
                        max_intensity[control.plate] = {}
                        med_intensity[control.plate] = {}

                    if channel not in intensity[control.plate].keys():
                        intensity[control.plate][channel] = []
                        max_intensity[control.plate][channel] = []
                        med_intensity[control.plate][channel] = []

                    intensity[control.plate][channel].append(np.sum(img[channel,]))
                    max_intensity[control.plate][channel].append(np.max(img[channel,]))
                    med_intensity[control.plate][channel].append(np.quantile(img[channel,][mask>0], 0.5))
        pb.close()
            
        self.nobjects=nobjects
        self.object_size=object_size
        self.intensity=intensity
        self.max_intensity=max_intensity
        self.med_intensity=med_intensity
                
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate per channel scaling factors based on previously calculated intensity_stats.tsv")
    parser.add_argument('-c','--controls', help='TSV file with <plate> <well> <channel> <name> describing control wells to use per plate', required=True)
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('-m','--mask_dir', help='Base dir to cell masks organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('--mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py cell. <pattern>.tiff", default="*_cell_mask_*_cp_masks.tiff")
    parser.add_argument('-o','--output', help='Output folder', default="./")
    parser.add_argument('-p','--plate', help='Plate(s) (subfolders) to process', nargs='+', default=None)
    #parser.add_argument('--plate_groups', help='File describing how plates are grouped, used for multicycle runs to ensure scaling is done per cycles. Use registration_manifest.tsv from tglow pipeline', default=None)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--whitelist', help='TSV file with "<plate>  <well>" on each row descrbing what to include.', default=None)
    parser.add_argument('--max_project', help="Calculate flatfields on max projections. Automatically activates --all_planes", action='store_true', default=False)
    #parser.add_argument('--q1', help='Quantile 1, the quantile in an image', default="q99.9999")
    #parser.add_argument('--q2', help='Quantile 2, the quantile over all images in input 0-100', default="95")
    #parser.add_argument('--pattern', help='File pattern of sumstats tsv', default="intensity_stats.tsv")
    #parser.add_argument('--scale_max', help='The max possible value of the input images, this will be the max value of the output', default=65535)
    args = parser.parse_args()

    scaler = ScalingCalculator(controls=args.controls,
                      input=args.input,
                      plates=args.plate,
                      masks=args.mask_dir,
                      output=args.output,
                      blacklist=args.blacklist,
                      mask_pattern=args.mask_pattern,
                      max_project=args.max_project)


    scaler.fetch_object_intensity()
    scaler.make_intensity_summary()
