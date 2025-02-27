#!/usr/bin/env python 

import os
import re
import tifffile
import string
import argparse
import random
import time
import numpy as np
import logging
import cv2 as cv
from basicpy import BaSiC
from matplotlib import pyplot as plt
from scipy import ndimage as ndi
from tglow.io.compound_image_provider import CompoundImageProvider
from tglow.utils.tglow_utils import float_to_16bit_unint
from skimage.filters import threshold_otsu, gaussian
from skimage.transform import downscale_local_mean, resize
from numpy.polynomial import polynomial as P

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


# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class FlatFieldTrainer():
    
    def __init__(self, path, output_dir, output_prefix,  channel,  nimg, merge_n, tune, fit_darkfield, all_planes=False, max_project=False, plates=None, fields=None, planes=None, threshold=False, autosegment=False, plot=False, blacklist=None, pseudoreplicates=0, blur=False, sigma=2):
        
        # Class to read images
        self.provider = CompoundImageProvider(path, nimg,  channel, blacklist, plates, fields, planes, pseudoreplicates, merge_n, max_project, all_planes)
        self.channel = channel
        
        self.nimg=nimg
        self.merge_n=merge_n
        self.all_planes=all_planes
        self.tune=tune
        self.fit_darkfield=fit_darkfield
        
        self.output_dir=output_dir
        self.output_prefix=output_prefix
    
        self.threshold=threshold
        self.autosegment=autosegment
        self.plot=plot
                
        self.blur=blur
        self.sigma=sigma
        
        # Output
        self.out = f"{self.output_dir}/{self.output_prefix}_ch{str(self.channel)}"
        
            
    def train_basicpy(self):
        
        # Read random images possibly as compound
        training_imgs = self.provider.fetch_training_images()
        
        # Apply a gaussian blur to the images prior to running bp
        if self.blur:
            i=0
            while i < len(training_imgs): 
                training_imgs[i] = gaussian(training_imgs[i], sigma=self.sigma, preserve_range=True)
                i+=1
                
        # Gather the weights
        if self.threshold:
            weights=[]
            k=0
            for img in training_imgs:
                threshold = threshold_otsu(img)
                log.debug(f"Threshold for img {k} = {threshold}")
                weights.append(img>threshold)
                k = k+1
            weights = np.array(weights)
        else:
            weights = None
    
        # Convert into 3d numpy array
        merged = np.array(training_imgs)
        log.info(f"Read training files into array of shape {str(merged.shape)}")

        # Init basicpy object with
        basic = BaSiC(get_darkfield=self.fit_darkfield,
                      autosegment=self.autosegment,
                      working_size=128)

        log.debug("Tuning parameters:")
        log.debug(basic.resize_params)
        log.debug(basic.working_size)

        # Optimze parameters
        if self.tune:
            log.info("Tuning model")

            init_params = {"smoothness_flatfield": 2}
            search_space = {"smoothness_flatfield": list(np.logspace(-0.5, 1.5, 15))}

            if self.fit_darkfield:
                init_params.update({
                    "smoothness_darkfield":  1e-3, 
                    "sparse_cost_darkfield":  1e-3
                })
                
                search_space.update({
                    "smoothness_darkfield":  list(np.logspace(-0.5, 1.5, 15)),
                    "sparse_cost_darkfield":  list(np.logspace(-0.5, 1.5, 15))            
                })
            
            basic.autotune(merged,
                           fitting_weight=weights,
                           #n_iter=200,
                           #histogram_qmin=0.05,
                           #histogram_qmax=0.95,
                           early_stop_torelance=1e-7,
                           search_space=search_space,
                           init_params=init_params,
                           early_stop_n_iter_no_change=25)
        
        # Fit parameters
        log.info("Fitting model")
        basic.fit(merged,
                fitting_weight=weights)
        
        if not os.path.exists(self.out):
            os.makedirs(self.out)
    
        basic.save_model(self.out, overwrite=True)

        if self.plot:
            self.plot_results(basic, merged)
    
    def train_mean_filter(self):
        
        # Read random images possibly as compound
        training_imgs = self.provider.fetch_training_images()
        
        sum_image=None
        sum_mask=None

        for final_img in training_imgs:
            final_img = final_img.astype(np.float32) 
            final_img = final_img / np.iinfo(np.uint16).max
            
            if sum_image is None:
                sum_image = np.zeros_like(final_img, dtype=np.float32)
                sum_mask = np.ones_like(final_img, dtype=np.float32)

                sum_image = sum_image + final_img      
                
                if self.threshold: 
                    thresh = threshold_otsu(final_img)
                    sum_mask += (final_img > thresh)

        if self.threshold: 
            avg = sum_image / sum_mask
        else:
            avg = sum_image / len(training_imgs)
        
        # Convert to unit8
        avg_norm = avg / np.max(avg)
        avg_norm = avg_norm * np.iinfo(np.uint8).max
        avg_norm = avg_norm.astype(np.uint8)
        
        blurred = cv.GaussianBlur(avg_norm, (1001,1001),500)
        
        # Normalize to mean to derive scaling factors
        flatfield = blurred / np.mean(flatfield)
        
        self.__create_output(flatfield,  np.array(training_imgs))
        
    def train_polynomial(self, degree):
    
        # Read random images possibly as compound
        training_imgs = self.provider.fetch_training_images()
        
        fit_sum = None
        i =0
        for img in training_imgs:
            log.info(f"Fitting poly to image {i}")
            cur_img = downscale_local_mean(img, (4, 4))
            
            if fit_sum is None:
                fit_sum = np.ones_like(cur_img)
            fit_sum += self.__fit_poly_img(cur_img, degree)
            i += 1
            
        # Calculate the average over the images
        fit_avg = fit_sum / len(training_imgs)
        final_fit = resize(fit_avg, training_imgs[0].shape, order=0)

        # Normalize flatfield to mean
        flatfield = final_fit / np.mean(final_fit)
        
        self.__create_output(flatfield, np.array(training_imgs))
    
    def __fit_poly_img(self, image, degree):
        Y, X = np.mgrid[:image.shape[0], :image.shape[1]]
        Z = image

        x = X.flatten()
        y = Y.flatten()
        z = Z.flatten()

        design_matrix = P.polyvander2d(x, y, [degree, degree])

        coeffs, residuals, rank, s = np.linalg.lstsq(design_matrix, z, rcond=None)
        coeffs = coeffs.reshape((degree + 1, degree + 1))

        z_fit = P.polyval2d(X, Y, coeffs)

        return(z_fit)

    
    def __create_output(self, flatfield, merged):
        # Hack, but save the results as a basicpy object to keep it compatible with the pipeline
        basic = BaSiC()
        basic.flatfield = flatfield
        basic.darkfield = np.zeros_like(flatfield)
        basic.baseline = np.mean(merged, axis=(1,2))
        
        if not os.path.exists(self.out):
            os.makedirs(self.out)
    
        basic.save_model(self.out, overwrite=True)
        
        if self.plot:
            self.plot_results(basic, merged)
        
        
    #def apply_transform(self, basic, images):
    #    im_float = images.astype(np.float32)
    #    output = (im_float - basic.darkfield[np.newaxis]) / basic.flatfield[np.newaxis]
    #    return(output)
        
    def plot_results(self, basic, merged):
        
        # Apply correction
        merged_corrected = basic.transform(merged)

        # Convert data back to original 16bit uint
        merged_corrected = float_to_16bit_unint(merged_corrected)
        
        # Plot results
        if self.plot:
            plot_basic_results(basic, self.out + "/flat_and_darkfield.png")
            plot_before_after_mp(np.max(merged, axis=0), np.max(merged_corrected, axis=0), filename=self.out + "/all_imgs_max_proj_pre_post.png")
            
            # Plot before and after for first 5 images, or as many as are available
            plot_max = 5
            if merged_corrected.shape[0] < plot_max:
                plot_max = merged_corrected.shape[0]
            
            i = 0
            while i < plot_max:
                plot_before_after(merged, merged_corrected, i, filename=self.out + f"/img{i}_pre_post.png")
                i += 1

if __name__ == "__main__":
    
    # CLI 
    parser = argparse.ArgumentParser(description="Train a basicpy model on raw HCI images orgnaized into <plate>/<row>/<col>/<field>.ome.tiff stacks with CZYX")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff')
    parser.add_argument('-o','--output', help='Output folder')
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--output_prefix', help='Output prefix name', default="flatfields")
    parser.add_argument('--no_tune', help="[BP only] Do not tune the basicpy model", action='store_true', default=False)
    parser.add_argument('--fit_darkfield', help="[BP only] Fit the darkfield component. For tglow not reccomended", action='store_true', default=False)
    parser.add_argument('--nimg', help="Number of random images to train on. Sampled with replacement", default=None, required=True)
    parser.add_argument('--merge_n', help="Number of images to combine into one --nimg times", default=1)
    parser.add_argument('--pseudoreplicates', help="Number of pseudoreplicates to generate. Pseudoreplicate is a compound from --merge_n images samples from --nimg read images", default=0)
    parser.add_argument('-p','--plate', help='Plate to process. Defaults to all detected plates.', nargs='+', default=None)
    parser.add_argument('-c','--channel', help="Channel number to correct", required=True)
    parser.add_argument('--fields', help='Fields to use. Defaults to use all fields.', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. Defaults to use all planes', nargs='+', default=None)
    #parser.add_argument('--gpu', help="Use the GPU", action='store_true', default=False)
    parser.add_argument('--all_planes', help="Instead of randomly picking one plane for a stack, use them all", action='store_true', default=False)
    parser.add_argument('--max_project', help="Calculate flatfields on max projections. Automatically activates --all_planes", action='store_true', default=False)
    parser.add_argument('--threshold', help="For mode MEAN, use otsu to only calculate mean in foreground regions. For BP Use Otsu threshold to set fitting_weight from basicpy to 1 for foreground and 0 for background. This conceptually the inverse of the autosegmentation", action='store_true', default=False)
    parser.add_argument('--autosegment', help="[BP only] Enable basicpy autosegment option", action='store_true', default=False)
    parser.add_argument('--blur', help="[BP only] Apply a gaussian blur prior to fitting model.", action='store_true', default=False)
    parser.add_argument('--sigma', help="[BP only] The sd of the gaussian kernel.", default=5)
    parser.add_argument('--plot', help="Plot basicpy results", action='store_true', default=False)
    parser.add_argument('--mode', help='The mode to run in BASICPY | MEAN | POLY', default="BASICPY")
    parser.add_argument('--degree', help="[POLY only] The degree for the polynomial", default=2)

    args = parser.parse_args()
    
    input=args.input
        
    # Set subfolders to process
    if args.plate == None:
        plates=[ name for name in os.listdir(input) if os.path.isdir(os.path.join(input, name)) ]
    else:
        plates=args.plate
    
    # If max projecting all_planes also needs to be true so the full stack is read
    # This is also forced in the constructor, but this is so its printed well here
    # TODO: make the printing a method
    if args.max_project:
        args.all_planes=True
    
    print("-----------------------------------------------------------")
    print("Input plates:\t" + str(plates))
    print("Input fields:\t" + str(args.fields))
    #print("Input planes:\t" + str(args.planes))    
    print("Input:\t\t" + input)
    print("Output:\t\t" + args.output)
    print("Output pre:\t" + args.output_prefix)
    print("channel:\t" + str(args.channel))
    print("n img:\t\t" + str(args.nimg))
    print("tune basicpy:\t" + str(not args.no_tune))
    print("merge n:\t" + str(args.merge_n))
    print("darkfied:\t" + str(args.fit_darkfield))
    print("max project:\t" + str(args.max_project))
    print("threshold:\t" + str(args.threshold))
    print("autosegment:\t" + str(args.autosegment))
    print("plot:\t" + str(args.plot))
    print("blacklist:\t" + str(args.blacklist))
    print("pseudoreplicates:\t" + str(args.pseudoreplicates))
    print("gaussian blur:\t" + str(args.blur))
    print("sigma:\t" + str(args.sigma))
    print("degree:\t" + str(args.degree))

    print("-----------------------------------------------------------")

    trainer = FlatFieldTrainer(path=input,
                            output_dir=args.output,
                            output_prefix=args.output_prefix,
                            channel=args.channel,
                            nimg=int(args.nimg),
                            tune=not args.no_tune,
                            merge_n=int(args.merge_n),
                            fit_darkfield=args.fit_darkfield,
                            all_planes=args.all_planes,
                            max_project=args.max_project,
                            plates=plates,
                            fields=args.fields,
                            threshold=args.threshold,
                            autosegment=args.autosegment,
                            plot=args.plot,
                            blacklist=args.blacklist,
                            pseudoreplicates=int(args.pseudoreplicates),
                            blur=args.blur,
                            sigma=float(args.sigma))
    
    if (args.mode == "BASICPY"):                 
        trainer.train_basicpy()
    elif (args.mode == "MEAN"):
        log.warning("MEAN method is very dodgy and I do not reccomend it, Use it at your own risk!")
        trainer.train_mean_filter()
    elif (args.mode == "POLY"):
        trainer.train_polynomial(int(args.degree))
    else:
        raise ValueError(f"{args.mode}, not a valid mode")






    
    
    
    