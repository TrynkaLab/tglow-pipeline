#!/usr/bin/env python 

from math import sqrt
from skimage import data, transform, io, exposure
from skimage.feature import blob_log
from skimage.color import rgb2gray
from skimage.filters import threshold_otsu, threshold_multiotsu
from skimage.registration import phase_cross_correlation
from scipy.spatial.distance import cdist
from scipy import ndimage
from tglow.io.tglow_io import AICSImageReader

import matplotlib.pyplot as plt
from scipy.ndimage import shift
from tifffile import imread, imwrite
from tglow.utils.tglow_plot import plot_grey_as_magma, plot_registration_imgs, composite_images, plot_blobs, plot_histogram
from tglow.utils.tglow_utils import float_to_16bit_unint
import numpy as np
import os

from pystackreg import StackReg
import pystackreg
import copy
import logging
import argparse

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#--------------------------------------------------------------
def filter_close_points(points, distance_threshold):
    # Calculate pairwise distances
    distances = cdist(points[:,0:2],points[:,0:2])
    
    # Create a mask for points that are close to each other
    # We use the upper triangle of the distance matrix to avoid self-comparisons
    #mask = np.triu(distances < distance_threshold, k=1)
    mask = distances < distance_threshold
    np.fill_diagonal(mask, False)
    
    # Find indices of points to keep (not close to any other point)
    indices_to_keep = ~np.any(mask, axis=1)
    
    # Apply the mask to filter the original array
    filtered_points = points[indices_to_keep]
    
    return filtered_points

#--------------------------------------------------------------
def find_extreme_values(data, percent=5):
    lower_percentile = np.percentile(data, percent)
    upper_percentile = np.percentile(data, 100 - percent)
    return [lower_percentile <= x <= upper_percentile for x in data]

#--------------------------------------------------------------
def filter_by_sd(data, sd_threshold=2):
    """
    Filter values within a specified standard deviation threshold and return their indexes.
    
    Parameters:
    - data: array-like, the input data
    - sd_threshold: float, the number of standard deviations to use as the threshold (default: 2)
    
    Returns:
    - numpy array of indexes where values are within the SD threshold
    """
    
    # Convert input to numpy array if it's not already
    data = np.array(data)
    
    # Calculate mean and standard deviation
    mean = np.mean(data)
    std = np.std(data)
    
    # Calculate lower and upper bounds
    lower_bound = mean - sd_threshold * std
    upper_bound = mean + sd_threshold * std
    
    # Find indexes where values are within the bounds
    indexes = np.where((data >= lower_bound) & (data <= upper_bound))[0]
    
    return indexes

#--------------------------------------------------------------
def fetch_blobcrops(raw_img, img_prefix, border_window, xy_window, bead_dist):

    # Remove blobs that have a in this window
    bead_window = round(bead_dist)
    
    # Define the window to consider blobs
    border_limits = {"ymin": border_window,
                    "ymax": raw_img.shape[1]-border_window,
                    "xmin": border_window,
                    "xmax": raw_img.shape[2]-border_window}

    # Max project
    raw_img_max = np.max(raw_img, axis=0)

    # Plot 
    #plot_grey_as_magma(raw_img_max, f"{outdir}/fields/f{img_prefix}_raw_img_max.png")

    norm_img_max = raw_img_max / np.max(raw_img_max)

    # Plot 
    #plot_grey_as_magma(norm_img_max, f"{outdir}/fields/f{img_prefix}_norm_img_max.png")

    # Identify blobs in 2d
    blobs_log = blob_log(norm_img_max,
                        min_sigma=1,
                        max_sigma=5,
                        num_sigma=10,
                        threshold=0.01,
                        exclude_border=True)

    # Calculate the radii of the blobs
    #blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)
    plot_blobs(norm_img_max, f"{outdir}/fields/f{img_prefix}_blobs_raw.png", blobs_log)
    np.savetxt(f'{outdir}/fields/f{img_prefix}_blobs_raw.tsv', blobs_log, delimiter='\t')

    log.info(f"Detected {blobs_log.shape[0]} beads")

    # Filter blobs closer together then the window
    blobs_filtered = filter_close_points(blobs_log, bead_window)

    # Filter XY
    keeps = (blobs_filtered[:,0] >= border_limits["ymin"]) & (blobs_filtered[:,0] <= border_limits["ymax"])
    keeps = keeps & (blobs_filtered[:, 1] >= border_limits["xmin"]) & (blobs_filtered[:,1] <= border_limits["xmax"])

    blobs_filtered = blobs_filtered[keeps,:]

    plot_blobs(norm_img_max, f"{outdir}/fields/f{img_prefix}_blobs_pos_filtered.png", blobs_filtered)
    np.savetxt(f'{outdir}/fields/f{img_prefix}_blobs_pos_filtered.tsv', blobs_filtered, delimiter='\t')

    log.info(f"Kept {blobs_filtered.shape[0]} beads after filtering for neighbours & position. {round((blobs_filtered.shape[0] / blobs_log.shape[0])*100, 2)}%")

    bead_crops = []
    for blob in blobs_filtered:
        y, x, r = blob
        
        # Y, if window is even add 1 so there is a center pixel
        ymin = round( y - xy_window)
        if (xy_window*2) % 2 == 0:
            ymax = round( y + xy_window) + 1
        else:
            ymax = round( y + xy_window)

        # X, if window is even add 1 so there is a center pixel
        xmin = round( x - xy_window)
        if (xy_window*2) % 2 == 0:
            xmax = round( x + xy_window) + 1
        else:
            xmax = round( x + xy_window)

        cur_crop=raw_img[:,ymin:ymax, xmin:xmax]
        cur_crop=cur_crop/np.max(cur_crop)
        bead_crops.append(cur_crop)
        
    return bead_crops

#--------------------------------------------------------------
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Create a PSF from a series of images")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff')
    parser.add_argument('-p','--plate', help='Plate to process. Defaults to all detected plates.', nargs='+', default=None)
    parser.add_argument('-w','--well', help='Wells to process (at least one). A01 [A02 B03 ...]', nargs='+', required=True)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    parser.add_argument('-c','--channel', help="Channel to distill psf from", required=True)
    parser.add_argument('-o','--output', help='Output folder')
    parser.add_argument('--suffix', help='Suffix to add to the output', default="")
    parser.add_argument('--window', help="xy distance to generate the PSF arround", default=16)
    parser.add_argument('--exclude', help="xy distance arround the edge of the field to exclude blobs", default=20)
    parser.add_argument('--max_distance', help="The distance to consider blobs to close together to use", default=20)
    args = parser.parse_args()

    #args.input = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/CRISPR_1/pipeline/results/images"
    #args.output = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/PSF_2/pipeline/results/testing"
    #args.plate = None
    #args.fields = None

    outdir = f"{args.output}"
    suffix = args.suffix
    xy_window = int(args.window)
    border_window = int(args.exclude)
    args.max_distance = int(args.max_distance)
    
    if not os.path.exists(f"{outdir}/fields"):
        os.makedirs(f"{outdir}/fields")

    if args.window >= args.exclude:
        raise RuntimeError("--window is larger then --exclude, this is not allowed as it means the crop will be over the edge of the image")

    #------------------------------------------------------------------
    # Read the bead images
    reader = AICSImageReader(args.input,
                              plates_filter=args.plate,
                              fields_filter=args.fields)
    
    # Get the indexed image queries
    avail_iqs = [x for xs in reader.images.values() for x in xs]    
    i = 1
    
    # Array to save the individual blobcrops
    bead_crops = []
    # Loop overal all images, and select the matching wells, read them, and fetch bead crops
    for iq in avail_iqs:
        for well in args.well:
            if (iq.get_well_id() == well):
                iq.channel = args.channel
                raw_img = reader.read_image(iq)   

                bead_crops.extend(fetch_blobcrops(raw_img, i, border_window, xy_window, args.max_distance))
                i+=1

    nblobs_raw = len(bead_crops)
    
    #------------------------------------------------------------------
    # Calculate the intensity of the beads to identify outliers
    bead_intensity = []
    bead_intensity_mean = []

    for cur_crop in bead_crops:
        bead_intensity.append(np.sum(cur_crop))
        
        # Calculate the mean intensity in the thresholded region
        thresh = threshold_otsu(cur_crop)
        tmp = copy.deepcopy(cur_crop)
        tmp[tmp < thresh] = 0
        
        bead_intensity_mean.append(np.sum(tmp) / np.sum(tmp > 0))

    # Plot histogram of bead intensities
    plot_histogram(bead_intensity, f"{outdir}/bead_integrated_intensity_whole_stack.png")
    plot_histogram(bead_intensity_mean, f"{outdir}/bead_mean_intensity_thresholded_region.png")

    beads_to_keep = np.array(filter_by_sd(bead_intensity_mean, 1.644854))

    # Filter beads which are SD outliers
    bead_intensity_mean = [bead_intensity_mean[i] for i in beads_to_keep]
    bead_intensity = [bead_intensity[i] for i in beads_to_keep]
    bead_crops = [bead_crops[i] for i in beads_to_keep]
    log.info(f"Kept {len(bead_intensity_mean)} beads after filtering for intensity. {round((len(bead_intensity_mean) / nblobs_raw)*100, 2)}%")

    #------------------------------------------------------------------
    # Registration
    beads_final = []
    ref = bead_crops[0]

    # Calculate the mean background signal using the outer 10% voxels in x and y
    background_mask = []
    background_mask.extend([ i for i in range(0, round(ref.shape[1]*0.1))])
    background_mask.extend([i for i in range((ref.shape[1] - round(ref.shape[1]*0.1)), ref.shape[1])])

    backgrounds=[]
    for img in bead_crops:
        backgrounds.append(np.percentile(img[:,background_mask, background_mask], 90))
    background = np.mean(backgrounds)

    errors=[]
    for img in bead_crops[1:len(bead_crops)]:
        reg_res = phase_cross_correlation(ref, img, upsample_factor=5)
        img_reg = shift(img, shift=reg_res[0])
        img_reg = img_reg - background
        img_reg[img_reg < 0] = 0
        img_reg = img_reg / np.max(img_reg)
        beads_final.append(img_reg)
        errors.append(reg_res[1])

    # Filter outliers based on RMS error
    beads_to_keep = np.array(filter_by_sd(errors, 1.644854))
    beads_final = [beads_final[i] for i in beads_to_keep]

    log.info(f"Kept {len(beads_final)} beads after filtering for registration error outliers. {round((len(beads_final) / nblobs_raw)*100, 2)}%")
    plot_histogram(errors, f"{outdir}/registration_rms_error.png")

    #------------------------------------------------------------------
    # Assemble final PSF
    beads_final = np.array(beads_final)
    # Its easy enough to rerun, so no need to save this
    #imwrite(f"{outdir}/registered_beads.tiff",  beads_final)

    # Aggregate mean of final PSF
    psf = np.mean(beads_final, axis=0)
    imwrite(f"{outdir}/psf_raw_means.tiff",  psf)

    # set values < 0 to 0
    psf_final = copy.deepcopy(psf)
    psf_final[psf_final < 0] = 0
    psf_final = psf_final / np.max(psf_final)

    # Align arround the center of mass
    center = np.round([psf.shape[0]/2, psf.shape[1]/2, psf.shape[1]/2])

    psf_cm = copy.deepcopy(psf_final)
    psf_cm[psf_cm < threshold_otsu(psf_cm)] = 0
    center_of_mass = ndimage.center_of_mass(psf_cm)

    psf_final = ndimage.shift(psf_final, center-center_of_mass)

    # Re normalize post shift
    psf_final = psf_final / np.max(psf_final)

    # Find the background region
    psf_final_background = np.max(psf_final[:,background_mask, background_mask])

    # Trim arround the edges to black borders that could be introduced during registration
    psf_final = psf_final[2:(psf_final.shape[0]-2),2:(psf_final.shape[1]-2),2:(psf_final.shape[2]-2)]

    psf_final_thresh = copy.deepcopy(psf_final)
    #thresh = threshold_multiotsu(psf_final, classes=2)[0]
    #psf_final_thresh[psf_final_thresh < thresh] = 0
    psf_final_thresh = psf_final_thresh - psf_final_background

    # Covert to uint16
    psf_final = psf_final / np.max(psf_final)
    psf_final = float_to_16bit_unint(psf_final * np.iinfo(np.uint16).max)

    psf_final_thresh = psf_final_thresh / np.max(psf_final_thresh)
    psf_final_thresh = float_to_16bit_unint(psf_final_thresh * np.iinfo(np.uint16).max)

    imwrite(f"{outdir}/psf_centered_bg_rm{suffix}.tiff",  psf_final_thresh)
    imwrite(f"{outdir}/psf_centered{suffix}.tiff", psf_final)

    # The below code saves a 5x and 6x downsampled version, this is no longer needed
    # #------------------------------------------------------------------
    # planes_to_keep=[]
    # planes_to_keep.append(round(psf_final_thresh.shape[0]/ 2))

    # for i in range(1, round(round(psf_final_thresh.shape[0]/6)/ 2)):
    #     planes_to_keep.append(round(psf_final_thresh.shape[0]/2) - (i*6))
    #     planes_to_keep.append(round(psf_final_thresh.shape[0]/2) + (i*6))

    # planes_to_keep.sort() 

    # psd_final_thresh_x = psf_final_thresh[planes_to_keep,:,]
    # psd_final_thresh_x = (psd_final_thresh_x / np.max(psd_final_thresh_x)) * np.iinfo(np.uint16).max
    # psd_final_thresh_x = float_to_16bit_unint(psd_final_thresh_x)

    # imwrite(f"{outdir}/psf_centered_bg_rm_6x_{suffix}.tiff",  psd_final_thresh_x)

    # #------------------------------------------------------------------
    # planes_to_keep=[]
    # planes_to_keep.append(round(psf_final_thresh.shape[0]/ 2))

    # for i in range(1, round(round(psf_final_thresh.shape[0]/5)/ 2)):
    #     planes_to_keep.append(round(psf_final_thresh.shape[0]/2) - (i*5))
    #     planes_to_keep.append(round(psf_final_thresh.shape[0]/2) + (i*5))

    # planes_to_keep.sort() 

    # psd_final_thresh_x = psf_final_thresh[planes_to_keep,:,]
    # psd_final_thresh_x = (psd_final_thresh_x / np.max(psd_final_thresh_x)) * np.iinfo(np.uint16).max
    # psd_final_thresh_x = float_to_16bit_unint(psd_final_thresh_x)

    # imwrite(f"{outdir}/psf_centered_bg_rm_5x_{suffix}.tiff",  psd_final_thresh_x)