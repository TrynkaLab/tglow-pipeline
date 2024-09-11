#!/usr/bin/env python 

# Adapted from Martin Prete's scripts
# Works with both Opera Phenix and Operetta indices
####################################################################
# Script invocation:
# python actions/parse_xml.py \
#    --input_file "/path/to/Images/Index.xml" \
#    --output_path "/path/to/output"
####################################################################
# Store Index.xml information in JSON format for easier parsing
# and re-use of the pate/well/images metadata.
####################################################################

import os
import re
import json
import logging
import argparse
from xml.etree import ElementTree as ET
from tglow.io.perkin_elmer_parser import PerkinElmerParser
#from tglow.io.image_query import ImageQuery

# Setup logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse Perkin Elmer index XML file to JSON file')
    parser.add_argument('--input_file', type=str, required=False, help='Path to the PE Index file')
    parser.add_argument('--output_path', type=str, required=True, help='Path to the output directory where to storer the JSON file')
    parser.add_argument('--to_manifest', required=False, action='store_true', default=False, help='Store the output as a <well> <plate> <index> tsv file')
    
    try:
       args = parser.parse_args()
    except:
       parser.print_help()
       exit(1)

    pe_data = PerkinElmerParser(args.input_file)
    
    if args.to_manifest:
        pe_data.write_manifest(args.output_path)
        pe_data.save(args.output_path)
        
        # Channel info
        info_file=open(f"{args.output_path}/acquisition_info.txt", 'w')
        
        info_file.write("#-----------------------------------------------------------------------------------\n")
        info_file.write("# Channel list:\n")
        if pe_data.channels is not None:
            
            i = 0
            info_file.write(f"index\tpe_id\tname\tnew_name\tresolution_zyx\texcitation\temission\texposure\n")
            for channel in pe_data.channels:
            
                res=channel['size']
                res = (len(pe_data.planes), res[1], res[0])
            
                info_file.write(f"{i}\t{channel['id']}\t{channel['name']}\tch{channel['id']} - {channel['name']}\t{res}\t{channel['main_excitation_wavelength']}\t{channel['main_emission_wavelength']}\t{channel['exposure_time']}\n")
                i+=1
        else:
            channel_names = None
            info_file.write("None detected\n")
            
        info_file.write("\n")
        
        # Objective info
        info_file.write("#-----------------------------------------------------------------------------------\n")
        info_file.write("# Objective info:\n")
        if pe_data.channels is not None:
            channel = pe_data.channels[0]
            info_file.write(f"acquisition_type: {channel['acquisition_type']}\n")
            info_file.write(f"binning_x: {channel['binning_x']}\n")
            info_file.write(f"binning_y: {channel['binning_y']}\n")
            info_file.write(f"objective_na: {channel['objective_na']}\n")
            info_file.write(f"objective_magnification: {channel['objective_magnification']}\n")
        else:
            channel_names = None
            info_file.write("None detected\n")
        
        info_file.write("\n")

        # ZYX pixel sizes
        resolution = pe_data.estimate_pixel_sizes()
        
        info_file.write("#-----------------------------------------------------------------------------------\n")
        info_file.write("# Pixel sizes: \n")
        if resolution is not None:
            info_file.write(f"ZYX: {resolution}\n")
            info_file.write(f"Estimated anisotropy: {resolution[0] / resolution[1]}\n")
        else:
            physical_pixel_sizes=None
            info_file.write("None detected\n")
                    
        info_file.flush()
        info_file.close() 
        
    else:
        pe_data.save(args.output_path)

