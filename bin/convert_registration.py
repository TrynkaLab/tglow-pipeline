#!/usr/bin/env python 

import os
import argparse
import numpy as np

import tglow.io.tglow_io as tglow_io
import pickle
import logging
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_utils import write_bin

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#------------------------------------------------------------------
# Main runner
class ConvertRegistration:
    
    # Constructor, and parse arguments
    def __init__(self, args):
          
        self.input=args.input
        self.output=args.output

        # Plates to process
        self.plate_reader=AICSImageReader(self.input, plates_filter=args.plate, fields_filter=args.fields, pattern="*.pickle")

        self.plates = [key for key in self.plate_reader.images.keys()]

        if args.fields is None:
            self.fields = self.plate_reader.fields[self.plates[0]]
            log.info(f"Autodetected fields {self.fields} based on first plate and image")
        else:
            self.fields=args.fields
            
    def fetch_registration(self, iq) -> dict:
        pickle_path = f"{self.input}/{iq.plate}/{ImageQuery.ID_TO_ROW[iq.row]}/{iq.col}/{iq.field}.pickle"
        
        if os.path.isfile(pickle_path):
            pickle_file = open(pickle_path, "rb")
            alignment_matrices = pickle.load(pickle_file)
            pickle_file.close()
        else:
            raise RuntimeError(f"Alignment matrix for {pickle_path} does not exist")

        return alignment_matrices  

    # Main loops
    def run(self):
        for plate in self.plates:
            log.info(f"Processing plate: {plate}")
            for cur_img in self.plate_reader.images[plate]:
                reg = self.fetch_registration(cur_img)
                for reg_plate in reg.keys():
                    
                    if reg_plate != "transform":
                        
                        outdir_final = f"{self.output}/{cur_img.plate}/{ImageQuery.ID_TO_ROW[cur_img.row]}/{cur_img.col}"
                        outfile = f"{outdir_final}/{cur_img.field}__{reg_plate}.bin"
                        
                        os.makedirs(outdir_final, exist_ok=True)
                        write_bin(reg[reg_plate], outfile)
       
    
    def printParams(self):
        
        log.info(f"Input:\t\t{self.input}")
        log.info(f"Output:\t\t{self.output}")    
        log.info(f"Input plates:\t{str(self.plates)}")
        log.info(f"Input fields:\t{str(self.fields)}")
 
        log.info(f"----------------------------")    




# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Convert registration output to a binary R compatible format")
    parser.add_argument('-i','--input', help='Base dir to  input')
    parser.add_argument('-o','--output', help='Output prefix')
    parser.add_argument('-p','--plate', help='Subfolder in dir to process. <plate1> | <plate1> <plate2> <plateN>', nargs='+', default=None)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    args = parser.parse_args()
  
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    
    runner=ConvertRegistration(args)
    runner.printParams()
    log.info("-----------------------------------------------------------")
    
    runner.run()

    #main(well, input, output, ref_channel, qry_channel, write_zstack, write_max_projection, basicpy_models, uint32)







