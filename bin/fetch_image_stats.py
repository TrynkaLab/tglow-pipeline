import numpy as np
import time
import logging
import argparse
import pandas as pd
from tglow.io.tglow_io import AICSImageReader, AICSImageWriter
from tglow.io.image_query import ImageQuery

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Simple utility to generate per well intensity_stats.tsv on a folder of existing images")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('-o','--output', help='Output folder', default="./")
    parser.add_argument('-p','--plates', help='Plates to process (at least one)', nargs='+', default=None)
    args = parser.parse_args()

    #args.input = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/CELL_DIVIDER/pipeline/results/images"
    #args.plates = None
    
    reader = AICSImageReader(args.input, plates_filter=args.plates)

    writer = AICSImageWriter(args.output)
    
    for plate in reader.index.keys():
        cur_plate = reader.index[plate]
        for row in cur_plate.keys():
            cur_row = cur_plate[row]
            for col in cur_row.keys():
                cur_col = cur_row[col]
                
                for iq in cur_col.values():
                    cur_stack = reader.read_stack(iq)
                    writer.write_stack(cur_stack, iq, stats_only=True)
                    
                writer.write_image_stats(iq)
