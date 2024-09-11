#!/usr/bin/env python 

import logging
import tifffile
import os
import argparse
import RedLionfishDeconv as rlf
import numpy as np
import pyopencl as cl

from clij2fft.richardson_lucy import richardson_lucy_nc, richardson_lucy
from tglow.io.tglow_io import AICSImageReader, BlacklistReader, AICSImageWriter
from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_utils import float_to_16bit_unint, float_to_32bit_unint, dict_to_str, float_to_16bit_unint_scaled

root_log = logging.getLogger()
root_log.setLevel(logging.INFO)

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class RunDecon:
    
    def __init__(self, args):
        
        self.input = args.input
        self.plate = args.plate
        self.well = args.well
        self.fields = args.fields
        
        self.output = args.output

        self.blacklist = args.blacklist
        self.niter = int(args.niter)
        self.regularization = float(args.regularization)
        
        self.mode = args.mode
        
        if self.mode.startswith("clij2"):
            try:
                platforms = cl.get_platforms()
                if len(platforms) <= 0:
                    raise RuntimeError("Could not find a valid open cl platform. Check your enviroment.")
                devices=platforms[0].get_devices()

                for device in devices:
                    log.info(f"Found open CL device: {device}")
                    log.info(f"Device has {device.get_info(cl.device_info.GLOBAL_MEM_SIZE)} mem available.")
            except:
                raise RuntimeError("Could not find a valid open cl platform. Check your enviroment.")
      
        self.max_project = args.max_project
        self.uint32 = args.uint32
        self.clip_max=int(args.clip_max)

        if self.blacklist is not None:
            bl_reader = BlacklistReader(self.blacklist)
            bl = bl_reader.read_blacklist()
        else:
            bl = None
        
        self.reader = AICSImageReader(self.input, plates_filter=args.plate, fields_filter=self.fields, blacklist=bl)
        self.writer = AICSImageWriter(self.output)
        
        if self.fields is None:
            self.fields = self.reader.fields[self.plate]
                    
        if args.channels is None:
            img = self.reader.get_img(self.reader.images[self.plate][0])
            self.channels = range(0, img.dims['C'][0])
        else:
            self.channels = [int(x) for x in args.channels]
        
        self.psf_string = args.psf
        
        # Read PSFs
        self.psfs = {}
        for val in args.psf:
            keypair = val.split("=")
            log.info(f"Reading PSF: {keypair}")
            self.psfs[int(keypair[0])] = tifffile.imread(keypair[1])
            log.info(f"Read PSF {keypair[0]} with shape: {self.psfs[int(keypair[0])].shape}")
            
        # Trim PSFs
        self.psf_stepsize=args.psf_subsample_z
        self.psf_planes=args.psf_crop_z

        for psf_key in self.psfs.keys():
            
            # Select every self.psf_stepsize planes out of the PSF
            if self.psf_stepsize is not None:
                
                self.psf_stepsize = int(self.psf_stepsize)
                psf_final = self.psfs[psf_key]
                planes_to_keep=[]
                planes_to_keep.append(round(psf_final.shape[0]/ 2))

                for i in range(1, round(round(psf_final.shape[0]/self.psf_stepsize)/ 2)):
                    planes_to_keep.append(round(psf_final.shape[0]/2) - (i*self.psf_stepsize))
                    planes_to_keep.append(round(psf_final.shape[0]/2) + (i*self.psf_stepsize))

                planes_to_keep.sort() 
                self.psfs[psf_key] = psf_final[planes_to_keep,:,]
                log.info(f"Selected every {self.psf_stepsize} planes from PSF {psf_key}. Shape: {self.psfs[psf_key].shape}")

            # Select the planes out of the PSF
            if self.psf_planes is not None:
                self.psf_planes = int(self.psf_planes)

                psf_final=self.psfs[psf_key]
                planes_to_keep=[]
                planes_to_keep.append(round(psf_final.shape[0]/ 2))

                for j in range(1, self.psf_planes):
                    planes_to_keep.append(round(psf_final.shape[0]/2) - j)
                    planes_to_keep.append(round(psf_final.shape[0]/2) + j)

                planes_to_keep.sort() 
                self.psfs[psf_key] = psf_final[planes_to_keep,:,]

                log.info(f"Selected {self.psf_planes} planes +- arround the middle of the stack for PSF {psf_key}. Shape {self.psfs[psf_key].shape}")
                
        #for psf_key in self.psfs.keys():
        #    tifffile.imwrite(f"psf_{psf_key}_final.tiff", self.psfs[psf_key])

    
    def run(self):
        
        for field in self.fields:
            self.run_decon(field)
    
    def run_decon(self, field):
        
        row, col = ImageQuery.well_id_to_index(self.well)
        iq = ImageQuery(self.plate, row, col, field)
        image = self.reader.read_image(iq)
        img = self.reader.get_img(iq)
        
        log.info(f"Read image of shape {image.shape}")

        decons = []
        for channel in self.channels:
            
            if (channel in self.psfs):  
                log.info(f"Deconvoluting field {field} channel {channel}")

                if self.mode == "clij2":
                    cur_decon = self.deconvolute_clij2(image, channel)
                elif self.mode == "clij2_nc":
                    cur_decon = self.deconvolute_clij2_nc(image, channel)
                elif self.mode.startsWith("rlf_"):
                    # This is already done by clij2
                    #pad_len = round(image.shape[1]/2)+1
                    #image = np.pad(image, ((0,0),(pad_len, pad_len), (0, 0), (0, 0)), mode='constant', constant_values=0)
                    #log.info(f"Padded image with zeroes in Z direction for final shape {image.shape}.")
                    cur_decon = self.deconvolute_rlf(image, channel)
                    # This is already done by clij2
                    #cur_decon = cur_decon[pad_len:-pad_len,:,:]
                else:
                    raise(f"{self.mode} is not a valid decon mode")
                
                log.debug(f"dtype: {cur_decon.dtype} max: {np.max(cur_decon)}")
            else:
                log.info(f"NOT deconvoluting field {field} channel {channel}")
                cur_decon = image[channel,:,:,:]
                
            if self.max_project:
                decons.append(np.max(cur_decon, axis=0, keepdims=True))
            else:
                decons.append(cur_decon)
        
        
        if self.uint32:        
            decon = float_to_32bit_unint(np.array(decons))
        else:
            decon = float_to_16bit_unint_scaled(np.array(decons), self.clip_max)
        
        log.debug(f"Final dtype: {decon.dtype} max: {np.max(decon)}")

        log.info(f"Writing stack of dims {decon.shape}")
        self.writer.write_stack(decon, iq, img.channel_names, img.physical_pixel_sizes)
        
        
    def deconvolute_clij2(self, image, channel):
    
        psf = self.psfs[channel]
        decon = richardson_lucy(image[channel,:,:,:],
                                psf,
                                numiterations=self.niter,
                                regularizationfactor=self.regularization)
        
        return(decon)
    
    def deconvolute_clij2_nc(self, image, channel):
    
        psf = self.psfs[channel]
        decon = richardson_lucy_nc(image[channel,:,:,:],
                                   psf,
                                   numiterations=self.niter,
                                   regularizationfactor=self.regularization)
        
        return(decon)
               
    def deconvolute_rlf(self, image, channel):
    
        psf = self.psfs[channel]
        decon = rlf.doRLDeconvolutionFromNpArrays(image[channel,:,:,:],
                                        psf,
                                        niter=self.niter,
                                        method= "gpu" if self.mode == "rlf_gpu" else "cpu")
        
        return(decon)
        
    
    def printParams(self):
        
        log.info(f"Input:\t\t{self.input}")
        log.info(f"Input plate:\t{str(self.plate)}")
        log.info(f"Input well:\t{str(self.well)}")
        log.info(f"Input fields:\t{str(self.fields)}")

        log.info(f"Output:\t\t{self.output}")
        log.info(f"psfs:\t\t{self.psf_string}")
        log.info(f"niter:\t\t{self.niter}")
        log.info(f"mode:\t\t{self.mode}")
        log.info(f"channels:\t\t{self.channels}")
        log.info(f"max project:\t\t{self.max_project}")
        log.info(f"Output 32bit:\t{str(self.uint32)}")    
        
     
# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Merge raw Phenix tiff files organized per well into per channel tiffs, one page per z-stack")
    parser.add_argument('-i','--input', help='Base dir to raw input')
    parser.add_argument('-o','--output', help='Output prefix')
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process')
    parser.add_argument('-w', '--well', help='Well ID to merge')
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    parser.add_argument('--psf', help='Point spread functions <zero index channel id>=/path/to/tiff', nargs='+')
    parser.add_argument('--niter', help='Number of RL iterations', default=10)
    parser.add_argument('--regularization', help='Regularization parameter. Only used for CLI2', default=0.0002)
    parser.add_argument('--mode', help='Mode of caluclating RL. clij2 | clij2_nc | rlf_cpu | rlf_gpu', default="clij2")
    parser.add_argument('--max_project', help="Output max projection over Z", action='store_true', default=False)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--clip_max', help='Clip values above this prior to 32>16 bit conversion. Values below this will be preserved, but scaled and rounded to 16 bit range', default=65535)
    parser.add_argument('--psf_crop_z', help='Number of planes arround the center of the PSF in z to crop. Defaults to all', default=None)
    parser.add_argument('--psf_subsample_z', help='Select every x planes out of the PSF. If psf has a Z spacing of 100nm and your data 500nm select 5. Defaults to all planes.', default=None)

    parser.add_argument('--channels', help='Channels to output, zero indexed. Only channels specified in --psf are deconveluted!', nargs='+', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint after applying basicpy model", action='store_true', default=False)
    args = parser.parse_args()
  
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    
    
    log.debug(f"FIELDS: {args.fields}")
    runner=RunDecon(args)
    runner.printParams()
    log.info("-----------------------------------------------------------")
    
    runner.run()