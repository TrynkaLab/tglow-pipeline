import numpy as np
import time
import logging
import argparse
import pandas as pd
import glob
import os

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate per channel scaling factors based on previously calculated intensity_stats.tsv")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('-o','--output', help='Output folder', default="./")
    parser.add_argument('--q1', help='Quantile 1, the quantile in an image', default="q99.9")
    parser.add_argument('--q2', help='Quantile 2, the quantile over all images in input 0-100', default="99")
    parser.add_argument('--pattern', help='File pattern of sumstats tsv', default="intensity_stats.tsv")
    parser.add_argument('--scale_max', help='The max possible value of the input images, this will be the max value of the output', default=65535)
    args = parser.parse_args()
    
    path = args.input
    
    args.q2 = float(args.q2)
    args.scale_max = float(args.scale_max)

    files = glob.glob(f"{path}/**/{args.pattern}", recursive=True)
    
    main_df = None
    
    for file in files:
        log.debug(f"Reading: {file}")
        cur_df = pd.read_csv(file, sep="\t")
        
        if main_df is not None:
            main_df = pd.concat((main_df, cur_df))
        else:
            main_df = cur_df

    channels = main_df["channel"].unique()

    scale_factors = []
    stats_max = []
    stats_mean = []
    stats_median = []
    
    quantiles = [0, 0.01, 0.1, 25, 50, 75, 99, 99.9, 100]

    for channel in channels: 
        cur_df = main_df[main_df["channel"] == channel]
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
        for plate in main_df["plate"].unique():
            scale_factors.append(f"{plate}_ch{channel}={scale_factor}")
    
    info_file=open(f"{args.output}/scaling_factors.txt", 'w')
    info_file.write(" ".join(scale_factors))
    info_file.flush()
    info_file.close()

    # Write the summary file
    header = ['type', 'channel', 'q0', 'q0.1', 'q1', 'q25', 'q50', 'q75', 'q99', 'q99.9', 'q100', 'mean']
    mean_file=open(f"{args.output}/intensity_summary.tsv", 'w')
    mean_file.write("\t".join(header))
    mean_file.write("\n")
    
    for ch in range(len(channels)):
        mean_file.write(f"mean_intensity\t{channels[ch]}\t")
        mean_file.write("\t".join([str(num) for num in stats_mean[ch]]))
        mean_file.write("\n")
        
        mean_file.write(f"median_intensity\t{channels[ch]}\t")
        mean_file.write("\t".join([str(num) for num in stats_median[ch]]))
        mean_file.write("\n")
    
        mean_file.write(f"max_intensity\t{channels[ch]}\t")
        mean_file.write("\t".join([str(num) for num in stats_max[ch]]))
        mean_file.write("\n")
    
    mean_file.flush()
    mean_file.close()