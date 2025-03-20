#!/usr/bin/env python 

import os
import argparse
import csv
import logging
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery


# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Main loop
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Create a manifest.tsv for a plate/row/col/fields.ome.tiff")
    parser.add_argument('-i','--input', help='Path to input root plate/row/col/fields.ome.tiff', default="./")
    parser.add_argument('-o','--output', help='Output prefix', default="./")
    parser.add_argument('-p','--plates', help='Subfolder in raw dir to process', nargs='+', default=None)
    parser.add_argument('--pattern', help="The pattern to index fields.", default="*.ome.tiff")
    args = parser.parse_args()
    
    # args.input = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/processed_images"
    # args.output = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/processed_images"
    # args.pattern = "*.ome.tiff"
    # args.plates=None
    
    reader = AICSImageReader(path=args.input, plates_filter=args.plates, pattern=args.pattern)
    
    for plate in reader.images.keys():
        cur_plate = reader.wells[plate]
        
        if not os.path.exists(f"{args.output}/{plate}"):
            os.makedirs(f"{args.output}/{plate}")
            
        outfile =  open(f"{args.output}/{plate}/manifest.tsv", 'w')
        
        writer = csv.writer(outfile, delimiter="\t")
        for well in cur_plate:
            row, col = ImageQuery.well_id_to_index(well)
            writer.writerow([well, ImageQuery.ID_TO_ROW[str(row)], col, plate, "none"])
            
        outfile.flush()
        outfile.close()