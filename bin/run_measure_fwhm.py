#!/usr/bin/env python

from inspect import stack
import numpy as np
import tifffile as tiff
import pandas as pd
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
import argparse
import logging
import matplotlib.pyplot as plt
import os
from tglow.utils.tglow_plot import plot_grey_as_magma, plot_registration_imgs, composite_images, plot_blobs, plot_histogram, plot_histogram_df
from skimage.registration import phase_cross_correlation
from scipy.ndimage import center_of_mass, shift
import copy
from skimage.filters import threshold_otsu
import tqdm
import random
from scipy.optimize import curve_fit
from scipy.ndimage import maximum_filter1d

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

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


def mad_scores(x):
    # x: iterable of numbers
    x = list(x)
    med = np.nanmedian(x)
    abs_devs = [abs(v - med) for v in x]
    mad = np.nanmedian(abs_devs)
    if mad == 0:
        # Avoid division by zero (all values equal or extremely tight)
        return [0.0 for _ in x]
    return [(v - med) / mad for v in x]

def mask_by_mad(series, mad_threshold=2.5):
    # indexes to keep within this column
    idx_keep = filter_by_mad(series.values, mad_threshold=mad_threshold)
    mask = np.zeros(len(series), dtype=bool)
    mask[idx_keep] = True
    return mask

def filter_by_mad(data, mad_threshold=2.5):
    """
    Filter values within a specified MAD threshold and return their indexes.
    
    Parameters:
    - data: array-like, the input data
    - mad_threshold: float, the number of MADs to use as the threshold (default: 2.5)
    
    Returns:
    - numpy array of indexes where values are within the MAD threshold
    """
    
    # Convert input to numpy array if it's not already
    data = np.array(data)
    data = np.array(mad_scores(data))
    
    # Calculate lower and upper bounds (typically use 2.5 * MAD for ~99% coverage)
    lower_bound = -mad_threshold
    upper_bound = mad_threshold
    
    # Find indexes where values are within the bounds
    indexes = np.where((data >= lower_bound) & (data <= upper_bound))[0]
    
    return indexes


def plot_bead_slices_gaussian(stack, filename, pixel_width=0):
    """
    Plot and save orthogonal slices + Gaussian-fitted profiles to specified filename.
    Z-axis is vertical (y-axis) for XZ and YZ plots. Viridis colormap.
    """
    # Use basename for titles
    base_name = os.path.splitext(os.path.basename(filename))[0]
    z0, y0, x0 = stack.shape[0]//2, stack.shape[1]//2, stack.shape[2]//2
    
    # === Orthogonal Slices ===
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # XY slice
    ax1.imshow(stack[z0], cmap='viridis', origin='lower')
    #ax1.plot(x0, y0, 'c*', markersize=15, markeredgewidth=2)
    ax1.set_title(f'{base_name}\nXY slice (Z={z0})')
    ax1.axis('tight')
    plt.colorbar(ax1.images[0], ax=ax1)
    
    # XZ slice - Z vertical (y-axis), X horizontal (x-axis)
    ax2.imshow(stack[:, y0, :], cmap='viridis', origin='lower', 
               extent=[0, stack.shape[2], 0, stack.shape[0]])
    #ax2.plot(x0, z0, 'c*', markersize=15, markeredgewidth=2)
    ax2.set_title(f'{base_name}\nXZ slice (Y={y0})')
    ax2.set_xlabel('X'); ax2.set_ylabel('Z')
    plt.colorbar(ax2.images[0], ax=ax2)
    
    # YZ slice - Z vertical (y-axis), Y horizontal (x-axis)
    ax3.imshow(stack[:, :, x0], cmap='viridis', origin='lower',
               extent=[0, stack.shape[1], 0, stack.shape[0]])
    #ax3.plot(y0, z0, 'c*', markersize=15, markeredgewidth=2)
    ax3.set_title(f'{base_name}\nYZ slice (X={x0})')
    ax3.set_xlabel('Y'); ax3.set_ylabel('Z')
    plt.colorbar(ax3.images[0], ax=ax3)
    
    plt.suptitle(f'PSF Orthogonal Slices - {base_name}', fontsize=16)
    plt.tight_layout()
    filename_slices = filename.replace('.png', '_slices') if '.png' in filename else f'{filename}_slices'
    #plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{filename_slices}.pdf", format='pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    # === Profiles with Gaussian Fits ===
    fig, (axz, axy, axx) = plt.subplots(1, 3, figsize=(15, 5))
    
    def extract_profile(stack, axis, pixel_width):
        z0, y0, x0 = stack.shape[0]//2, stack.shape[1]//2, stack.shape[2]//2
        half_width = pixel_width
        if axis == 0:
            slice_y = min(half_width, y0, stack.shape[1]-y0-1)
            slice_x = min(half_width, x0, stack.shape[2]-x0-1)
            return np.mean(stack[:, y0-slice_y:y0+slice_y+1, x0-slice_x:x0+slice_x+1], axis=(1,2))
        elif axis == 1:
            slice_z = min(half_width, z0, stack.shape[0]-z0-1)
            slice_x = min(half_width, x0, stack.shape[2]-x0-1)
            return np.mean(stack[z0-slice_z:z0+slice_z+1, :, x0-slice_x:x0+slice_x+1], axis=(0,2))
        else:
            slice_z = min(half_width, z0, stack.shape[0]-z0-1)
            slice_y = min(half_width, y0, stack.shape[1]-y0-1)
            return np.mean(stack[z0-slice_z:z0+slice_z+1, y0-slice_y:y0+slice_y+1, :], axis=(0,1))
    
    def fit_gaussian_profile(profile):
        """Fit Gaussian and return params, fitted curve, or None if fails."""
        if len(profile) < 5 or np.max(profile) <= 0:
            return None
        
        x = np.arange(len(profile))
        max_idx = np.argmax(profile)
        height_guess = profile[max_idx]
        center_guess = x[max_idx]
        sigma_guess = 2.0
        
        try:
            bounds = ([0, x[0], 0, -np.inf], [np.inf, x[-1], np.inf, np.inf])
            popt, _ = curve_fit(gaussian, x, profile, 
                              p0=[height_guess, center_guess, sigma_guess, 0],
                              bounds=bounds, maxfev=1000)
            fitted = gaussian(x, *popt)
            fwhm = 2.35482 * popt[2]  # FWHM = 2*sqrt(2*ln2)*sigma
            return popt, fitted, fwhm
        except:
            return None
    
    for axis, ax, title in [(0, axz, 'Z-profile'), (1, axy, 'Y-profile'), (2, axx, 'X-profile')]:
        profile = extract_profile(stack, axis, pixel_width)
        pos = np.arange(len(profile))
        
        # Plot data
        ax.plot(pos, profile, 'r-o', linewidth=2, markersize=4, label='Data')
        ax.axvline(len(profile)//2, color='c', linestyle='--', alpha=0.7, label='Center')
        
        # Gaussian fit
        fit_result = fit_gaussian_profile(profile)
        if fit_result is not None:
            popt, fitted, fwhm = fit_result
            ax.plot(pos, fitted, 'b-', linewidth=2, label=f'Gaussian fit\nFWHM={fwhm:.2f}px')
            # Mark half-max on fit
            half_max = fitted[np.argmax(fitted)] / 2
            ax.axhline(half_max, color='g', linestyle=':', alpha=0.8)
        else:
            ax.text(0.5, 0.5, 'Fit failed', transform=ax.transAxes, ha='center', va='center')
        
        ax.set_title(f'{base_name}: {title}')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Bead profile + Gaussian Fits (pixel_width={pixel_width}) - {base_name}', fontsize=14)
    plt.tight_layout()
    profiles_fn = filename.replace('.png', '_profiles') if '.png' in filename else f'{filename}_profiles'
    #plt.savefig(f"{profiles_fn}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{profiles_fn}.pdf", format='pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    print(f"Plots saved:\n- Slices: {filename}\n- Profiles+Gaussian: {profiles_fn}")

def get_bead_crops_com(blobs_filtered, raw_img, xy_window):
    """
    Extract bead crops from raw image based on filtered blob coordinates,
    center them using center of mass, and return the list of centered bead crops
    along with the mean background signal.
        
    :param blobs_filtered: ZYX positions of where to crop the beads
    :param raw_img: the image array with the beads, zyx
    :param xy_window: window size in pixels to crop around each bead center
    """
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
        #cur_crop=cur_crop/np.max(cur_crop)
        bead_crops.append(cur_crop)
         
    log.info(f"Found {len(bead_crops)} raw beads in images")
    
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

    log.info("Centering bead crops on the center of mass")
    
    progress = tqdm.tqdm(total=len(bead_crops)-1, desc="Centering beads", unit="bead")
    for img in bead_crops[1:len(bead_crops)]:
        progress.update(1)
        
        img_cm = copy.deepcopy(img)
        img_cm[img_cm < threshold_otsu(img_cm)] = 0
        
        # Calculate center of mass (z, y, x order)
        com = center_of_mass(img_cm)
                
        # Target center (middle of stack)
        center = np.round([img.shape[0]/2, img.shape[1]/2, img.shape[2]/2])

        # Shift stack to center (order=3 for cubic interpolation)
        centered = shift(img, center-com)
        beads_final.append(centered)
    
    return beads_final, background

def fwhm(stack, axis, pixel_width=0, mode='emperical'):
    """
    Calculate FWHM along specified axis from stack center.
    Assumes PSF peak is centered in middle of ZYX stack.
    
    Parameters:
    - stack: np.ndarray shape (Z, Y, X) 
    - axis: 0=Z, 1=Y, 2=X
    - pixel_width: int, pixels added to each side of center
    
    Returns: FWHM in pixels
    """
    
    axes_coords = [np.arange(s) for s in np.array(stack.shape)]
    
    z0, y0, x0 = stack.shape[0]//2, stack.shape[1]//2, stack.shape[2]//2
    half_width = pixel_width
    
    # Extract line profile with averaging window
    if axis == 0:  # Z-axis
        slice_y = min(half_width, y0, stack.shape[1]-y0-1)
        slice_x = min(half_width, x0, stack.shape[2]-x0-1)
        profile = np.mean(stack[:, 
                               y0-slice_y:y0+slice_y+1, 
                               x0-slice_x:x0+slice_x+1], 
                         axis=(1,2))
        
    elif axis == 1:  # Y-axis  
        slice_z = min(half_width, z0, stack.shape[0]-z0-1)
        slice_x = min(half_width, x0, stack.shape[2]-x0-1)
        profile = np.mean(stack[z0-slice_z:z0+slice_z+1, :, x0-slice_x:x0+slice_x+1], 
                         axis=(0,2))
        
    else:  # X-axis
        slice_z = min(half_width, z0, stack.shape[0]-z0-1)
        slice_y = min(half_width, y0, stack.shape[1]-y0-1)
        profile = np.mean(stack[z0-slice_z:z0+slice_z+1, 
                               y0-slice_y:y0+slice_y+1, :], 
                         axis=(0,1))
    
    
    if mode == 'emperical':
        # FWHM calculation with subpixel interpolation
        half_max = np.max(profile) / 2
        peak_idx = np.argmax(profile)
        
        # Left crossing
        left_idxs = np.where(profile[:peak_idx] >= half_max)[0]
        if len(left_idxs) == 0:
            return np.nan
        left_idx = left_idxs[0]
        if left_idx > 0:
            left_pos = left_idx - (half_max - profile[left_idx-1]) / (profile[left_idx] - profile[left_idx-1])
        else:
            left_pos = 0
            
        # Right crossing  
        right_idxs = np.where(profile[peak_idx:] >= half_max)[0]
        if len(right_idxs) == 0:
            return np.nan
        right_idx = peak_idx + right_idxs[-1]
        if right_idx < len(profile) - 1:
            right_pos = right_idx + (half_max - profile[right_idx]) / (profile[right_idx+1] - profile[right_idx])
        else:
            right_pos = len(profile) - 1    
        
        return right_pos - left_pos
    elif mode == 'gaussian':
        return fit_profile(profile, axes_coords[axis], init_sigma=2)
    else:
        raise ValueError("Invalid mode. Choose 'emperical' or 'gaussian'.") 
    
def gaussian(x, height, mu, sigma, offset=0):
    """1D Gaussian function: offset + height * exp(-(x-mu)^2 / (2*sigma^2))"""
    return offset + height * np.exp(-(x - mu)**2 / (2 * sigma**2))

def fit_profile(profile, x, init_sigma=2):
    """Fit Gaussian and return FWHM."""
    if len(profile) < 5 or np.max(profile) <= 0:
        return np.nan
    
    max_idx = np.argmax(profile)
    init_height = profile[max_idx]
    init_mu = x[max_idx]
    init_offset = 0

    try:
        # Bounds prevent nonsensical fits
        bounds = ([0, x[0], 0, -np.inf], [np.inf, x[-1], np.inf, np.inf])
        popt, _ = curve_fit(gaussian, x, profile, 
                            p0=[init_height, init_mu, init_sigma, init_offset],
                            bounds=bounds, maxfev=1000)
        sigma = popt[2]
        return 2.35482 * sigma, popt  # FWHM = 2 * sqrt(2*ln(2)) * sigma
    except:
        return np.nan, np.nan



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Evaluate the FWHM of beads in images.")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff')
    parser.add_argument('-p','--plate', help='Plate to process. Defaults to all detected plates.', nargs='+', default=None)
    parser.add_argument('-w','--well', help='Wells to process (at least one). A01 [A02 B03 ...]', nargs='+', required=True)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    parser.add_argument('-c','--channel', help="Channel to evaluate", required=True)
    parser.add_argument('-o','--output', help='Output folder')
    parser.add_argument('--window', help="xy distance to generate the PSF arround", default=16)
    parser.add_argument('--blobs', help="Path to the fields folder from run_beads_to_psf.py output", default=16)
    parser.add_argument('--pixel_size', help="Pixel size in microns", nargs='+', default=None)
    args = parser.parse_args()


    #------------------------------------------------------------------
    args.channel = int(args.channel)
    if args.pixel_size is not None:
        args.pixel_size = [float(x) for x in args.pixel_size]
        if len(args.pixel_size) != 3:
            raise ValueError("Pixel size must have 3 values: z y x in microns")
    
    
    #------------------------------------------------------------------
    os.makedirs(args.output, exist_ok=True)
    
    # Read the bead images
    reader = AICSImageReader(args.input,
                              plates_filter=args.plate,
                              fields_filter=args.fields)
    
    # Get the indexed image queries
    avail_iqs = [x for xs in reader.images.values() for x in xs]    
    i = 1
    
    # Array to save the individual blobcrops
    final_beads = []
    
    # Loop overal all images, and select the matching wells, read them, and fetch bead crops
    # Beads are registered against the first bead found
    for iq in avail_iqs:
        for well in args.well:
            if (iq.get_well_id() == well):
                
                #iq = ImageQuery(iq.plate, 5, 2, 1)
                iq.channel = args.channel
                log.info(f"Reading image: {iq.to_string()}")
                raw_img = reader.read_image(iq)   
                
                # Read previous filtered blobs used for making PSF, these should be valid 
                # and non-overlapping.
                blobpath = f"{args.blobs}/f{iq.field}_blobs_pos_filtered.tsv"
                blobs_df = pd.read_csv(blobpath, sep="\t", header=None, names=['y','x','r'])
                blobs_filtered = blobs_df[['y','x','r']].to_numpy()
                    
                beads, background = get_bead_crops_com(blobs_filtered, raw_img, args.window)
                final_beads.extend(beads)
                i += 1
    
    #--------------------------------------------
    # Calculate the intensity of the beads to identify outliers
    bead_intensity_mean = []

    for cur_crop in final_beads:   
        bead_intensity_mean.append(np.sum(cur_crop) / cur_crop.size)     

    # Plot histogram of bead intensities
    plot_histogram(bead_intensity_mean, f"{args.output}/stack_mean_intensity.png")
    
    # Filter beads by MAD
    beads_to_keep = filter_by_mad(bead_intensity_mean, mad_threshold=3)
    
    log.info(f"Removing {len(final_beads) - len(beads_to_keep)} beads based on intensity filtering (>3 MADs)")
    
    final_beads = [final_beads[i] for i in beads_to_keep]
    bead_intensity_mean = [bead_intensity_mean[i] for i in beads_to_keep]
    
    # Plot histogram of filtered bead intensities
    plot_histogram(bead_intensity_mean, f"{args.output}/stack_mean_intensity_filtered.png")
    
    #--------------------------------------------
    # Calculate FWHM for all beads using the emperical & gaussian approach
    log.info("Calculating FWHM for all beads")
    result = pd.DataFrame(columns=['fwhm_emp_z_px', 'fwhm_emp_y_px', 'fwhm_emp_x_px','fwhm_gau_z_px', 'fwhm_gau_y_px', 'fwhm_gau_x_px'])

    progress = tqdm.tqdm(total=len(final_beads)-1, desc="fwhm", unit="bead")
    for bead in final_beads:
        progress.update(1)
        
        row = len(result)
        result.loc[row, ['fwhm_emp_z_px', 'fwhm_emp_y_px', 'fwhm_emp_x_px']] = [fwhm(bead, axis=0), fwhm(bead, axis=1), fwhm(bead, axis=2)]
        
        for axis, col_prefix in [(0, 'fwhm_gau_z_px'), (1, 'fwhm_gau_y_px'), (2, 'fwhm_gau_x_px')]:
            fwhm_val, popt = fwhm(bead, axis=axis, mode='gaussian')
            result.loc[row, [col_prefix]] = fwhm_val

    log.info("Unfiltered mean FWHM values (pixels):")
    print(result[['fwhm_emp_z_px', 'fwhm_emp_y_px', 'fwhm_emp_x_px', 'fwhm_gau_z_px', 'fwhm_gau_y_px', 'fwhm_gau_x_px']].mean())
    
    # Calculate the resolution in microns if pixel size is provided
    if args.pixel_size is not None:
        result['fwhm_emp_z_micron'] = result['fwhm_emp_z_px'] * args.pixel_size[0]
        result['fwhm_emp_y_micron'] = result['fwhm_emp_y_px'] * args.pixel_size[1]
        result['fwhm_emp_x_micron'] = result['fwhm_emp_x_px'] * args.pixel_size[2]
        result['fwhm_gau_z_micron'] = result['fwhm_gau_z_px'] * args.pixel_size[0]
        result['fwhm_gau_y_micron'] = result['fwhm_gau_y_px'] * args.pixel_size[1]
        result['fwhm_gau_x_micron'] = result['fwhm_gau_x_px'] * args.pixel_size[2]
        log.info("Unfiltered mean FWHM values (microns):")
        print(result[['fwhm_emp_z_micron', 'fwhm_emp_y_micron', 'fwhm_emp_x_micron', 'fwhm_gau_z_micron', 'fwhm_gau_y_micron', 'fwhm_gau_x_micron']].mean())
    
    # Ensure numeric values
    for col in result.columns:
        result[col] = pd.to_numeric(result[col], errors='coerce')  # NaN for bad values

    #------------------------------------------------------------------
    # Remove outliers > 3 MAD deviations from the median for gaussian XYZ FWHM columns
    cols_to_check = ['fwhm_gau_z_px', 'fwhm_gau_y_px', 'fwhm_gau_x_px']

    # Calulate the combined mask over the 3 columns
    mask = np.ones(len(result), dtype=bool)
    for c in cols_to_check:
        mask &= mask_by_mad(result[c], mad_threshold=3)

    n_removed = (~mask).sum()
    log.info(f"Removing {n_removed} outlier rows (>3 MADs) from FWHM results")

    # Filter the beads and results
    result_filtered = result[mask]
    final_beads_filtered = [final_beads[i] for i in range(len(final_beads)) if mask[i]]

    result.to_csv(f"{args.output}/fwhm_values.tsv", sep="\t", index=False)
    result_filtered.to_csv(f"{args.output}/fwhm_values_filtered.tsv", sep="\t", index=False)

    if args.pixel_size is not None:
        plot_histogram_df(result_filtered[['fwhm_gau_z_micron', 'fwhm_gau_y_micron', 'fwhm_gau_x_micron']],
                        filename=f"{args.output}/fwhm_gau_histogram_micron",
                        title_prefix="FWHM Gaussian",
                        bins=50,
                        xlabel=f"FWHM microns - channel {args.channel}",
                        ylabel="Number of beads")
    else:
        plot_histogram_df(result_filtered[['fwhm_gau_z_px', 'fwhm_gau_y_px', 'fwhm_gau_x_px']],
                        filename=f"{args.output}/fwhm_gau_histogram_pixels",
                        title_prefix="FWHM Gaussian",
                        bins=50,
                        xlabel=f"FWHM pixels - channel {args.channel}",
                        ylabel="Number of beads")
        
    #--------------------------------------------
    # Plot 5 random bead profiles
    sample = random.sample(range(0, len(final_beads_filtered)), min(5, len(final_beads_filtered)))
    for i in sample:
        plot_bead_slices_gaussian(final_beads_filtered[i], pixel_width=0, filename=f"{args.output}/bead_{i}")

    
    