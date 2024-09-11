import os
import re
import tifffile
import string
import argparse
import json
import numpy as np
from basicpy import BaSiC
from pathlib import Path
from skimage import transform, io, exposure, color, util, filters, measure
from matplotlib import pyplot as plt
from pystackreg import StackReg

import pickle
import copy
import logging
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_plot import composite_images, plot_registration_imgs

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#------------------------------------------------------------------
# Main runner
class Registration:
    
    # Constructor, and parse arguments
    def __init__(self, args):
          
        self.input=args.input
        self.output=args.output

        # Plates to process
        self.ref_plate=args.plate
        if len(self.ref_plate) != 1:
            raise RuntimeError("Must only supply one plate for registration as the reference plate")
        self.ref_plate=args.plate[0]
        
        # Channels for alignment and plates to merge
        self.plates_merge=args.plate_merge   

        if self.plates_merge is not None:
            # Reader for plates
            self.plate_reader=AICSImageReader(self.input, plates_filter=[self.ref_plate] + self.plates_merge, fields_filter=args.fields)
        
            if (args.ref_channel is None) | (args.qry_channel is None):
                raise RuntimeError("Need to supply reference and query channels for registration")
            
            self.ref_channel=int(args.ref_channel[0])
            
            if not args.ref_channel_eval is None:
                self.ref_channel_eval = int(args.ref_channel_eval[0])
            
            #self.qry_channel=int(args.qry_channel[0])
            self.qry_channel={}
            self.qry_channel_eval={}
             
            i=0
            for cur_plate in self.plates_merge:
                self.qry_channel[cur_plate] = int(args.qry_channel[i])
                
                if not args.qry_channel_eval is None:
                    self.qry_channel_eval[cur_plate] = int(args.qry_channel_eval[i])
                i+=1      
        else:
            raise RuntimeError("Need to supply --plates_merge for running registration")
            
        
        # Registration options
        self.plot=args.plot
        self.save_reg=True
        self.transform=StackReg.TRANSLATION
        self.eval_merge=args.eval_merge
        
        if self.eval_merge:
            if (args.ref_channel_eval is None) | (args.qry_channel_eval is None):
                raise RuntimeError("Need to supply reference and query evaluation channels for registration evaluation")      
            
        # Check number of fields matches
        self.fields = self.plate_reader.fields[self.ref_plate]
        
        for plate_merge in self.plates_merge:
            
            if len(self.plate_reader.fields[plate_merge]) != len(self.fields):
                raise RuntimeError("Number of fields in merge plate must match number of fields in reference plate")     
        
    # Main loops
    def run(self, well):
    
        row, col = ImageQuery.well_id_to_index(well)
        
         # Construct the output dir
        reg_dir=f"{self.output}/{self.ref_plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}/"

        alignment_matrices={}
        
        if not os.path.exists(reg_dir):
            os.makedirs(reg_dir)
            log.info(f"Folder created: {reg_dir}")

        if self.eval_merge:
            log.info(f"Calculating pearson correlation between ref and qry for channels {self.ref_channel_eval} and {self.qry_channel_eval}")
            outfile = f"{reg_dir}/registration_eval_pearson_correlations.tsv"
            file = open(outfile, "w")
            file.write("well\trow\tcol\tfield\tref_plate\tqry_plate\tref_ch_eval\tqry_ch_eval\tref_ch\tqry_ch\tr_before\tp_before\tr_after\tt_after\n")
            file.close()
            
        log.info(f"Running alignment between plates on max projections")
                
        for field in self.fields:
            cur_iq = ImageQuery(self.ref_plate, row, col, field)
            
            alignment_matrices = self.align_field(cur_iq, reg_dir)
            if self.save_reg:
                alignment_matrices['transform'] = self.transform
                pickle.dump(alignment_matrices, open(f"{reg_dir}/{field}.pickle", "wb")) 
                
            if self.eval_merge:
                self.eval_field(cur_iq, alignment_matrices, outfile)
        
        
    # Perform registration between planes
    def align_field(self, iq, reg_dir):
            
        ref_iq = ImageQuery(self.ref_plate, iq.row, iq.col, iq.field, self.ref_channel)
        
        #log.debug(ref_iq.to_string())
        ref_stack = self.plate_reader.read_image(ref_iq)
    
        # Always max project, as we don't want different alignments between planes
        ref_stack = np.max(ref_stack, axis=0)
        
        alignment_matrices={}
        i=0
        # Loop over plates to merge onto ref
        for plate_merge in self.plates_merge:
            log.info(f"Processing field {iq.field} plate {plate_merge}")
            if i==0:
                alignment_matrices[plate_merge] = {} 
                i+=1
                
            qry_iq = ImageQuery(plate_merge, iq.row, iq.col, iq.field, self.qry_channel[plate_merge])
            
            qry_stack = self.plate_reader.read_image(qry_iq)
            
            # Always max project, as we don't want different alignments between planes
            qry_stack = np.max(qry_stack, axis=0)
            
            # Calculate the offsets
            sr = StackReg(self.transform)
            align_mat = sr.register(ref_stack, qry_stack)
                            
            # Save the offsets
            alignment_matrices[plate_merge] = align_mat

            # Plot 
            if self.plot:
                # Apply the offesets
                qry_stack_reg = sr.transform(qry_stack, tmat=align_mat)
                before_reg = composite_images([ref_stack, qry_stack])
                after_reg = composite_images([ref_stack, qry_stack_reg])
                
                outfile = f"{reg_dir}/{iq.field}_{plate_merge}_refch{self.ref_channel}_qrych{self.qry_channel[plate_merge]}.png"
                plot_registration_imgs(before_reg, after_reg, filename=outfile)
        
        return alignment_matrices

                
    def eval_field(self, iq, alignment_mats, outfile):
        
        for plate_merge in self.plates_merge:
            
            # Calculate the correlations between the evaluation channels before and after registering
            a = np.max(self.plate_reader.read_image(ImageQuery(self.ref_plate, iq.row, iq.col, iq.field, self.ref_channel_eval)), axis=0)
            c = np.max(self.plate_reader.read_image(ImageQuery(plate_merge, iq.row, iq.col, iq.field, self.qry_channel_eval[plate_merge])), axis=0)  
            
            sr = StackReg(self.transform)
            b = sr.transform(c, tmat=alignment_mats[plate_merge])
            
            apcc, apval = measure.pearson_corr_coeff(a, c)
            log.info(f"Before Pearson R: {apcc:0.3g}, p-val: {apval:0.3g}")

            bpcc, bpval = measure.pearson_corr_coeff(a, b)
            log.info(f"After Pearson R {bpcc:0.3g}, p-val: {bpval:0.3g}")
            
            line=f"{iq.get_well_id()}\t{iq.row}\t{iq.col}\t{iq.field}\t{self.ref_plate}\t{plate_merge}\t{self.ref_channel_eval}\t{self.qry_channel_eval[plate_merge]}\t{self.ref_channel}\t{self.qry_channel[plate_merge]}\t{apcc}\t{apval}\t{bpcc}\t{bpval}\n"
            file = open(outfile, 'a')
            file.write(line)
            file.close()

    # Load a previously generated alignment matrix in {plate: <matrix>} dict
    def load_field(iq, reg_dir):
        
        pickle_path = f"{reg_dir}/{iq.plate}/{ImageQuery.ID_TO_ROW[iq.row]}/{iq.col}/{iq.field}.pickle"
        alignment_matrices = pickle.load(pickle_path)

        return alignment_matrices
    
    def printParams(self):
        
        log.info(f"Input plates:\t{str(self.ref_plate)}")
        log.info(f"Input fields:\t{str(self.fields)}")
        #log.info(f"Input planes:\t{str(self.planes)}")   
        #log.info(f"Input channel:\t{str(self.channels)}")     
        log.info(f"Input:\t\t{self.input}")
        log.info(f"Output:\t\t{self.output}")    

        log.info(f"----------------------------")    
        log.info(f"Merge plates:\t{str(self.plates_merge)}")     
        #log.info(f"Merge channels:\t{str(self.channels_merge)}")  
        log.info(f"Ref channel:\t{str(self.ref_channel)}")  
        log.info(f"Qry channel:\t{str(self.qry_channel)}")
        log.info(f"Ref channel eval:\t{str(self.ref_channel_eval)}")  
        log.info(f"Qry channel eval:\t{str(self.qry_channel_eval)}")  
        log.info(f"Transform:\t{str(self.transform)}")        
        log.info(f"Save registr:\t{str(self.save_reg)}")     
        log.info(f"Plot:\t\t{str(self.plot)}")     


# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Merge raw Phenix tiff files organized per well into per channel tiffs, one page per z-stack")
    parser.add_argument('-w', '--well', help='Well ID to merge')
    parser.add_argument('-i','--input', help='Base dir to raw input')
    parser.add_argument('-o','--output', help='Output prefix')
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+')
    parser.add_argument('-r','--ref_channel', help='Reference channel for registration (nucleus)', nargs=1, default=None)
    parser.add_argument('-q','--qry_channel', help='Query channel for registration (nucleus)', nargs="+", default=None)
    parser.add_argument('-m','--plate_merge', help='Plates to align using pystackreg (e.g. multiple cycles of the same plate). Output will be added as additional channels.', nargs='+', default=None)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. <plane #1> | [<plane #1> <plane #2> <plane #n>]', nargs='+', default=None)
    parser.add_argument('--plot', help="[OPTIONAL] Plot the before and after registration", action='store_true', default=False)
    parser.add_argument('--eval_merge', help="[OPTIONAL] Calculate correlation between two channels before and after", action='store_true', default=False)
    parser.add_argument('--ref_channel_eval', help='[OPTIONAL] Reference channel for registration evaluation (nucleus)', nargs=1, default=None)
    parser.add_argument('--qry_channel_eval', help='[OPTIONAL] Query channel for registration evaluation (nucleus)', nargs='+', default=None)
    args = parser.parse_args()
  
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    
    runner = Registration(args)
    log.info(f"Well:\t{args.well}")
    runner.printParams()
    log.info("-----------------------------------------------------------")
    
    runner.run(args.well)

    #main(well, input, output, ref_channel, qry_channel, write_zstack, write_max_projection, basicpy_models, uint32)







