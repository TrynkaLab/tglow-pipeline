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
from basicpy import BaSiC
from matplotlib import pyplot as plt
from scipy import ndimage as ndi
from tglow.io.tglow_io import AICSImageReader, BlacklistReader
from tglow.utils.tglow_utils import float_to_16bit_unint
from skimage.filters import threshold_otsu, gaussian

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

class BasicpyTrainer():
    
    def __init__(self, path, output_dir, output_prefix,  channel,  nimg, merge_n, tune, fit_darkfield, all_planes=False, max_project=False, plates=None, fields=None, planes=None, threshold=False, autosegment=False, plot=False, blacklist=None, pseudoreplicates=0, blur=False, sigma=2):
        
        self.path=path
        self.channel=channel
        self.plates=plates
        self.fields=fields
        self.planes=planes
        
        self.nimg=nimg
        self.merge_n=merge_n
        self.all_planes=all_planes
        self.tune=tune
        self.fit_darkfield=fit_darkfield
        
        self.output_dir=output_dir
        self.output_prefix=output_prefix
        
        self.max_project=max_project
        
        if self.max_project:
            self.all_planes=True
            
        self.threshold=threshold
        
        self.autosegment=autosegment
        
        self.plot=plot
        
        self.blacklist = blacklist
        
        self.pseudoreplicates = pseudoreplicates
        
        self.blur=blur
        
        self.sigma=sigma
            
    def train(self):
        
        # Read the placklist
        if self.blacklist is not None:
            bl_reader = BlacklistReader(self.blacklist)
            bl = bl_reader.read_blacklist()
        else:
            bl = None
            
        # Init the image reader
        ac_reader = AICSImageReader(self.path, self.plates, self.fields, bl)
        training_imgs = []
        i=0

        start_time = time.time()
        total_imgs = 0
        
        if self.pseudoreplicates > 0:
            log.info(f"Pseudoreplicating {self.pseudoreplicates} compound images from {self.nimg} read images, with each compound image consisting of {self.merge_n} images")
            sample_n = 1
        else:
            sample_n = self.merge_n
        
        while i < self.nimg:
            i=i+1
            # Select subset of images
            
            all_imgs = [x for v in ac_reader.images.values() for x in v]
            
            merge_files = random.choices(all_imgs, k=sample_n)
            stack_dims = ac_reader.get_img(all_imgs[0]).dims
            
            #log.info(f"Reading {len(merge_files)} randomly selected files")

            # Loop over the raw tiffs and read them as a list of 2d numpy arrays
            j = 0
            training_imgs_tmp = []
            while j < len(merge_files): 
                #print(merge_files[j])
                #training_imgs_tmp.append(tifffile.imread(merge_files[j]))
                q = merge_files[j]
                q.channel = self.channel
                
                if self.all_planes:
                    img = ac_reader.read_image(q)
                    if self.max_project:
                        training_imgs_tmp.append(np.max(img, axis=0))
                    else:
                        training_imgs_tmp.append(list(img))
                else:
                    if self.planes is not None:
                        rplane = random.sample(self.planes, k=1)[0]
                    else:
                        rplane = random.sample(range(0, stack_dims["Z"][0]), k=1)[0]
                    q.plane = rplane
                    img = ac_reader.read_image(q)
                    training_imgs_tmp.append(img)
                log.debug(f"Read stack of {img.shape} [{total_imgs +1}/{len(merge_files) * self.nimg}]")
                j += 1
                total_imgs +=1
                
            if self.merge_n > 1:
                #log.info(f"Merging  {str(self.merge_n)} random images for training")
                training_imgs.append(np.max(np.array(training_imgs_tmp), axis=0))
            else:
                # This gave issues as it puts an array of arrays, giving the wrong shape
                #training_imgs.append(training_imgs_tmp)
                training_imgs=training_imgs_tmp
              
        log.info(f"Reading took { round(((time.time() - start_time)/60), 2)} minutes")  
        
        
        # Generate pseudoreplicates
        if self.pseudoreplicates > 0:
            log.info(f"Creating {self.pseudoreplicates} pseudoreplicates from {self.nimg} images")
            training_imgs_tmp=[]
            
            i = 0
            while i < self.pseudoreplicates:
                merge_files = random.choices(training_imgs, k=self.merge_n)
                training_imgs_tmp.append(np.max(np.array(merge_files), axis=0))
                i+=1
            
            training_imgs = training_imgs_tmp
            
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

        # Set background pixels to zero
        #if self.threshold:
        #    log.info("Calculating Otsu threshold prior to running basicpy")
        #    thresh = threshold_otsu(merged)
        #    #log.info(f"Setting pixels with value < {thresh} to zero")
        #    #merged[merged < thresh] = 0   
        #    weights=merged > thresh
        #else:
        #    weights=None
        
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
        
        # Save output
        out = f"{self.output_dir}/{self.output_prefix}_ch{str(self.channel)}"
        
        if not os.path.exists(out):
            os.makedirs(out)
            
        basic.save_model(out, overwrite=True)
        
        # Apply correction
        merged_corrected = basic.transform(merged)
        
        # Convert data back to original 16bit uint
        merged_corrected = float_to_16bit_unint(merged_corrected)
        
        # Plot results
        
        if self.plot:
            plot_basic_results(basic, out + "/flat_and_darkfield.png")
            plot_before_after_mp(np.max(merged, axis=0), np.max(merged_corrected, axis=0), filename=out + "/all_imgs_max_proj_pre_post.png")
            
            # Plot before and after for first 5 images, or as many as are available
            plot_max = 5
            if merged_corrected.shape[0] < plot_max:
                plot_max = merged_corrected.shape[0]
            
            i = 0
            while i < plot_max:
                plot_before_after(merged, merged_corrected, i, filename=out + f"/img{i}_pre_post.png")
                i += 1


if __name__ == "__main__":
    
    # CLI 
    parser = argparse.ArgumentParser(description="Train a basicpy model on raw HCI images orgnaized into <plate>/<row>/<col>/<field>.ome.tiff stacks with CZYX")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff')
    parser.add_argument('-o','--output', help='Output folder')
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--output_prefix', help='Output prefix name', default="basicpy")
    parser.add_argument('--no_tune', help="Do not tune the basicpy model", action='store_true', default=False)
    parser.add_argument('--fit_darkfield', help="Fit the darkfield component. For tglow not reccomended", action='store_true', default=False)
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
    parser.add_argument('--threshold', help="Use Otsu threshold to set fitting_weight from basicpy to 1 for foreground and 0 for background. This conceptually the inverse of the autosegmentation", action='store_true', default=False)
    parser.add_argument('--autosegment', help="Enable basicpy autosegment option", action='store_true', default=False)
    parser.add_argument('--plot', help="Plot basicpy results", action='store_true', default=False)
    parser.add_argument('--blur', help="Apply a gaussian blur prior to fitting model.", action='store_true', default=False)
    parser.add_argument('--sigma', help="The sd of the gaussian kernel.", default=5)

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

    print("-----------------------------------------------------------")

    log.warning("[DEPRECATED] This runner has been deprecated in favor of run_flatfield_estimation.py")

    trainer = BasicpyTrainer(path=input,
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
                            
    trainer.train()
    log.warning("[DEPRECATED] This runner has been deprecated in favor of run_flatfield_estimation.py")






    
    
    
    