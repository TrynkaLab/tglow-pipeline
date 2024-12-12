#!/usr/bin/env python 

import numpy as np
import time
import logging
import argparse
import pandas as pd
import glob
import os
from tglow.io.tglow_io import BlacklistReader

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate per channel scaling factors based on previously calculated intensity_stats.tsv")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('-o','--output', help='Output folder', default="./")
    parser.add_argument('-p','--plate', help='Subfolders to process', nargs='+', default=None)
    parser.add_argument('--plate_groups', help='File describing how plates are grouped, used for multicycle runs to ensure scaling is done per cycles. Use registration_manifest.tsv from tglow pipeline', default=None)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--q1', help='Quantile 1, the quantile in an image', default="q99.9999")
    parser.add_argument('--q2', help='Quantile 2, the quantile over all images in input 0-100', default="95")
    parser.add_argument('--pattern', help='File pattern of sumstats tsv', default="intensity_stats.tsv")
    parser.add_argument('--scale_max', help='The max possible value of the input images, this will be the max value of the output', default=65535)
    args = parser.parse_args()
    
    path = args.input
    
    args.q2 = float(args.q2)
    args.scale_max = float(args.scale_max)

    # Build list of files
    files = glob.glob(f"{path}/**/{args.pattern}", recursive=True)
    
    log.info(f"Indexed {len(files)} files")
    # Read the blacklist
    if args.blacklist is not None:
        bl_reader = BlacklistReader(args.blacklist)
        bl = bl_reader.read_blacklist_as_prc()
        log.info(f"Read blacklist with {len(bl)} patterns")
        log.info(f"Blacklist consists of patterns: {bl}")

        #files = [file for file in files if not reg]  
        files = [file for file in files if not any(pattern in file for pattern in bl)]
        log.info(f"Filtered using blacklist to {len(files)} files")

    if args.plate is not None:
        log.info(f"Platelist consists of patterns: {args.plate}")
        files = [file for file in files if any(pattern in file for pattern in args.plate)]
        log.info(f"Filtered using platelist to {len(files)} files")


    main_df = None

    indices = []
    scale_factors = []
    stats_max = []
    stats_mean = []
    stats_median = []
    quantiles = [0, 0.01, 0.1, 25, 50, 75, 95, 99, 99.9, 99.99, 99.999, 99.9999, 100]
    
    plate_groups = []
    
    # Define the plate groups (one cycle of imaging = a plate group)
    # This ensures channels are scaled witin a cycle and not accross cycles.
    if args.plate_groups != None:
        df = pd.read_csv(args.plate_groups, sep="\t")
        
        print(df.head())
                
        # Extract the first column (ref plate)
        plate_groups.append(df["reference_plate"].tolist())

        # Extract the third column (query plates)
        qry_plates = df["query_plates"].tolist()

        # Split the 3rd col on comma, and add each subsequent col
        # as a plategroup
        i = 0
        for qry in qry_plates:
            cur_qry = qry.split(",") 
            
            for j in range(0, len(cur_qry)):
                if i == 0:
                    plate_groups.append([])
                    
                plate_groups[j+1].append(cur_qry[j])
        
            i += 1
    
    else:
        plate_groups = [[main_df["plate"].unique()]]
    
    log.info("Considering plate groups:")
    log.info(plate_groups)
    
    
    log.info("Reading intensity files")
        # Read per well intensity CSV files
    for file in files:
        #log.debug(f"Reading: {file}")
        cur_df = pd.read_csv(file, sep="\t")
        
        if main_df is not None:
            main_df = pd.concat((main_df, cur_df))
        else:
            main_df = cur_df


    i = 0
    for plate_group in plate_groups:
        # Subset to the plates in the group to normalize
        plate_df = main_df[main_df["plate"].isin(plate_group)]

        # Find available channels
        channels = plate_df["channel"].unique()

        for channel in channels: 
            
            indices.append(f"group{i}_ch{channel}")
        
            # Subset to the current channel
            cur_df = plate_df[plate_df["channel"] == channel]
        
            # Calculate the scale factor
            scale_factor = np.percentile(cur_df[args.q1], args.q2) / args.scale_max
            
            # Max
            cur_stats = np.percentile(cur_df["q100"], quantiles).tolist()
            cur_stats.append(np.mean(cur_df["q100"]))
            stats_max.append(cur_stats)
            
            # Mean
            cur_stats = np.percentile(cur_df["mean"], quantiles).tolist()
            cur_stats.append(np.mean(cur_df["mean"]))
            stats_mean.append(cur_stats)
            
            # Median
            cur_stats = np.percentile(cur_df["q50"], quantiles).tolist()
            cur_stats.append(np.mean(cur_df["q50"]))
            stats_median.append(cur_stats)
            
            # Scale factors
            for plate in plate_df["plate"].unique():
                scale_factors.append(f"{plate}_ch{channel}={scale_factor}")
        i += 1
    
    # Write scaling factors
    info_file=open(f"{args.output}/scaling_factors.txt", 'w')
    info_file.write(" ".join(scale_factors))
    info_file.flush()
    info_file.close()

    # Write the summary file
    header = ['type', 'channel', 'q0', 'q0.1', 'q1', 'q25', 'q50', 'q75', 'q95', 'q99', 'q99.9', 'q99.99', 'q99.999', 'q99.9999', 'q100', 'mean']
    mean_file=open(f"{args.output}/intensity_summary.tsv", 'w')
    mean_file.write("\t".join(header))
    mean_file.write("\n")
    
    for ch in range(len(indices)):
        mean_file.write(f"mean_intensity\t{indices[ch]}\t")
        mean_file.write("\t".join([str(num) for num in stats_mean[ch]]))
        mean_file.write("\n")
        
        mean_file.write(f"median_intensity\t{indices[ch]}\t")
        mean_file.write("\t".join([str(num) for num in stats_median[ch]]))
        mean_file.write("\n")
    
        mean_file.write(f"max_intensity\t{indices[ch]}\t")
        mean_file.write("\t".join([str(num) for num in stats_max[ch]]))
        mean_file.write("\n")
    
    mean_file.flush()
    mean_file.close()
    
    log.info("Done")