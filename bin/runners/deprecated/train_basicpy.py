import os
import re
import tifffile
import string
import argparse
import random
import numpy as np
from basicpy import BaSiC
from matplotlib import pyplot as plt
from scipy import ndimage as ndi

# Plot results from basicpy fit
def plot_basic_results(basic, filename):
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    im = axes[0].imshow(basic.flatfield)
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Flatfield")
    im = axes[1].imshow(basic.darkfield)
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Darkfield")
    axes[2].plot(basic.baseline)
    axes[2].set_xlabel("Frame")
    axes[2].set_ylabel("Baseline")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()

# Plot before and after of an imageset for a given index in the arrays
def plot_before_after(images, images_transformed, i, filename):
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow(images[i])
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Original")
    im = axes[1].imshow(images_transformed[i])
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Corrected")
    fig.suptitle(f"frame {i}")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()

# Plot before and after for two pre-picked images, usefull for max projections
def plot_before_after_mp(images,images_transformed, filename):
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow(images)
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Original")
    im = axes[1].imshow(images_transformed)
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Corrected")
    fig.suptitle(f"Max projected")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()

# Build dict linking letters to number
def build_well_index(upper=True):
    if (upper):
        chars = string.ascii_uppercase
    else:
        chars = string.ascii_lowercase
    nums = list(range(1, len(chars)))
    idx = dict(zip(chars,nums))
    return(idx)

# Function that converts wells 
def make_prefix(well, idx, prepend_zero=True):
    match = re.search("([a-zA-Z])(\\d+)", well)
    if match != None:
        row = match[1]
        col = int(match[2])
        row_idx = idx[row]
        if (prepend_zero):
            prefix = "r" + str(row_idx).zfill(2) + "c" + str(col).zfill(2)
        else:
            prefix = "r" + str(row_idx) + "c" + str(row)
        return prefix
    else:
        return None

# Convert a numpy 32bit float to a 16 bit int, clip values to zero and 65535
def float_to_16bit_unint(matrix):

    # Convert to 16 bit to keep consistency with output
    # Round to nearest int            
    matrix = np.rint(matrix)
    
    # Set negatives to zero
    matrix[matrix < 0] = 0
    
    # Clip values at the 16 bit max
    matrix[matrix > 65535] = 65535
    
    matrix = matrix.astype(np.uint16)
    
    return matrix
         
# Main loop
def train_basicpy(input_dir_raw, output_prefix, output_dir, channel, nimg=100, tune=True, merge_n=1, blur=False, blur_sigma=1, get_darkfield=False):
    
    # List with available files
    training_files=[]
    
    # Loop over possible subfolders
    for plate in PLATES:
        
        print("[INFO] Building index of images for: " + plate)
        cur_path_root = input_dir_raw + "/" + plate + "/"
        
        # Select the available wells
        wells = os.listdir(cur_path_root)
        wells = [well for well in wells if re.search("[A-Z]\d\d", well) != None]
        
        #print(wells)
        print(PLANES)
        print(FIELDS)

        # Build index of all available files
        for well in wells:
            cur_path = cur_path_root + well + "/"
            well_files = os.listdir(cur_path)
            prefix=make_prefix(well, COLLUMN_INDEX) 

            for field in FIELDS:
                for plane in PLANES:
                    search_pattern = prefix + FIELD_PREFIX + str(field).zfill(2) + PLANE_PREFIX + str(plane).zfill(2) + "-" + CHANN_PREFIX + str(channel) + "sk1fk1fl1.tiff"

                    cur_files = [cur_path + file for file in well_files if re.search(search_pattern, file) != None]
                    #print(len(cur_files))
                    if (len(cur_files) >= 1):
                        training_files.append(cur_files)
                    
    print(f"[INFO] Found {str(len(training_files))} training files")
    
    training_imgs = []
    i=0

    while i < merge_n:
        i=i+1
        # Select subset of images
        merge_files = random.sample(training_files, nimg)
        
        print(f"[INFO] Reading {len(merge_files)} randomly selected files")
        #print(merge_files)
        # Loop over the raw tiffs and read them as 2d numpy arrays
        j = 0
        training_imgs_tmp = []
        while j < len(merge_files): 
            #print(merge_files[j])
            training_imgs_tmp.append(tifffile.imread(merge_files[j]))
            j=j+1
            
        if merge_n > 1:
            print("[INFO] Merging " + str(nimg) + " random images for training")
            training_imgs.append(np.max(np.array(training_imgs_tmp), axis=0))
        else:
            # This gave issues as it puts an array of arrays, giving the wrong shape
            #training_imgs.append(training_imgs_tmp)
            training_imgs=training_imgs_tmp
            
    # Convert into 3d numpy array
    merged_orig = np.array(training_imgs)
    print(f"[INFO] Read training files into array of shape {str(merged_orig.shape)}")

    if blur:
        print("[INFO] Applying gaussian blur to images prior to fitting model")
        merged = ndi.gaussian_filter(merged_orig, sigma=(0, blur_sigma, blur_sigma))
    else:
        merged=merged_orig

    # Init basicpy object with autosegmentation (needed for 3d)
    basic = BaSiC(get_darkfield=get_darkfield, smoothness_flatfield=1, autosegment=False)

    # Optimze parameters
    if tune:
        basic.autotune(merged)
    
    # Fit parameters
    basic.fit(merged)
    
    # Save output
    out = output_dir + "/" + output_prefix + "_ch" + str(channel)
    
    if not os.path.exists(out):
        os.makedirs(out)
        
    basic.save_model(out, overwrite=True)
    
    # Apply correction
    merged_corrected = basic.transform(merged_orig)
    
    # Convert data back to original 16bit uint
    merged_corrected = float_to_16bit_unint(merged_corrected)
    
    # Plot results
    plot_basic_results(basic, out + "/flat_and_darkfield.png")
    plot_before_after_mp(np.max(merged_orig, axis=0), np.max(merged_corrected, axis=0), filename=out + "/all_imgs_max_proj_pre_post.png")
    plot_before_after(merged_orig, merged_corrected, 0, filename=out + "/img0_pre_post.png")
    plot_before_after(merged_orig, merged_corrected, 1, filename=out + "/img1_pre_post.png")

    if blur:
        merged_corrected_blur = basic.transform(merged)
        plot_before_after_mp(np.max(merged, axis=0), np.max(merged_corrected_blur, axis=0), filename=out + "/all_imgs_max_proj_pre_post_blur.png")


# Main loop, needs a lot of cleanup, too many globals
if __name__ == "__main__":
    
    # Hardcoded globals for now    
    # Dict to convert well letters to numbers
    COLLUMN_INDEX=build_well_index()
    FIELD_PREFIX="f"
    PLANE_PREFIX="p"
    CHANN_PREFIX="ch"
    
    # CLI 
    parser = argparse.ArgumentParser(description="Merge raw Phenix tiff files organized per well into per channel tiffs, one page per z-stack")
    parser.add_argument('-i','--input', help='Base dir to raw input organized <plate>/<well>/<images>.tiff')
    parser.add_argument('-o','--output', help='Output folder')
    parser.add_argument('--output_prefix', help='Output prefix name', default="basicpy")
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+')
    parser.add_argument('-c','--channel', help="Channel number to correct")
    parser.add_argument('--no_tune', help="Do not tune the basicpy model", action='store_true', default=False)
    parser.add_argument('--fit_darkfield', help="Fit the darkfield component. For tglow not reccomended", action='store_true', default=False)
    parser.add_argument('--nimg', help="Number of random images to train on. If --merge_n is specified this is the number of images that are max projected into one.", default=100)
    parser.add_argument('--merge_n', help="Number of times to max project random selection --nimg number of images. If > 1 basicpy is run on this number of max projected images.", default=1)
    parser.add_argument('--blur', help="Blur the images with a guassian filter prior to fitting", action='store_true', default=False)
    parser.add_argument('--blur_sigma', help="SD of the gaussian filter if --blur", default=1)
    parser.add_argument('--fields', help='Fields to use', nargs='+', default=[1,2,3,4,5,6,7,8,9])
    parser.add_argument('--planes', help='Z planes to use', nargs='+', default=[1,2,3,4,5])
    args = parser.parse_args()
    
    input=args.input
    output_dir=args.output
    output_prefix=args.output_prefix
    
    # Set subfolders to process
    if args.plate == None:
        PLATES=[ name for name in os.listdir(input) if os.path.isdir(os.path.join(input, name)) ]
    else:
        PLATES=args.plate
    
    channel=int(args.channel)
    nimg=int(args.nimg)
    tune=not args.no_tune
    merge_n=int(args.merge_n)
    fit_darkfield=args.fit_darkfield
    blur=args.blur
    blur_sigma=float(args.blur_sigma)

    FIELDS=args.fields
    PLANES=args.planes

    print("-----------------------------------------------------------")
    print("Input plates:\t" + str(PLATES))
    print("Input fields:\t" + str(FIELDS))
    print("Input planes:\t" + str(PLANES))    
    print("Input:\t\t" + input)
    print("Output:\t\t" + output_dir)
    print("Output pre:\t" + output_prefix)
    print("channel:\t" + str(channel))
    print("n img:\t\t" + str(nimg))
    print("tune basicpy:\t" + str(tune))
    print("merge n:\t" + str(merge_n))
    print("darkfied:\t" + str(fit_darkfield))
    print("blur:\t\t" + str(blur))
    print("blur sigma:\t" + str(blur_sigma))
    print("-----------------------------------------------------------")

    train_basicpy(input, output_prefix, output_dir, channel, nimg, tune, merge_n, blur, blur_sigma, fit_darkfield)
    





