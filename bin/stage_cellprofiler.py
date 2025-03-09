#!/usr/bin/env python 

import os
import tifffile
import argparse
import numpy as np
from basicpy import BaSiC

import logging
from tglow.io.image_query import ImageQuery
from tglow.io.processed_image_provider import ProcessedImageProvider

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#------------------------------------------------------------------
# Main runner
class MergeAndAlign:
    
    # Constructor, and parse arguments
    def __init__(self, args):
        self.input=args.input
        self.output=args.output
        self.write_zstack=not args.no_zstack
        self.write_max_projection_onefile=args.max_project_onefile
        self.write_max_projection=args.max_project
        self.uint32=args.uint32
        self.fields=args.fields
        self.plates=args.plate
        
        if (self.write_max_projection):
            self.write_zstack=False
            log.info("Writing max projection as individual files, not saving z-stacks")
        
        self.provider = ProcessedImageProvider(path=args.input,
                                               blacklist=args.blacklist,
                                               plate=args.plate,
                                               plate_merge=args.plate_merge,
                                               registration_dir=args.registration_dir,
                                               flatfields=args.flatfields,
                                               scaling_factors=args.scaling_factors,
                                               mask_channels=args.mask_channels,
                                               mask_dir=args.mask_dir,
                                               mask_pattern=args.mask_pattern,
                                               uint32=args.uint32)

    # Main loops
    def run(self, well):

        self.provider.write_channel_index(self.output)
        row, col = ImageQuery.well_id_to_index(well)

        # Loop over possible plates
        for plate in self.plates:
                
            row, col = ImageQuery.well_id_to_index(well)
            out_dir_final= f"{self.output}/{plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}"
            
            if not os.path.exists(out_dir_final):
                os.makedirs(out_dir_final)
                log.info(f"Folder created: {out_dir_final}")
                    
            # Assume all plates have the same fields
            if self.fields is None:
                fields = self.provider.plate_reader.fields[plate]
            else:
                fields = self.fields
                
            for field in fields:
 
                iq = ImageQuery(plate, row, col, field)
                stack = self.provider.fetch_image(iq)
                
                #-------------------------------------------------
                # Save the images per channel
                # Array to store the channels for the max projection
                max_projection = []
                
                for channel in range(0, stack.shape[0]):
                    log.info(f"{plate} well: {well} channel: {str(channel)} field: {field}: Saving.")

                    # Final output filename for channel
                    cur_out = f"{out_dir_final}/{field}_{plate}_{well}_ch{str(channel)}.tiff"      
                    merged = stack[channel,:,:,:]

                    # Write out the final tiff file, saves as ZYX
                    if self.write_zstack:
                        if self.uint32:
                            tifffile.imwrite(cur_out, merged, shape=merged.shape, photometric='MINISBLACK', metadata={'axes': 'ZYX'})
                        else:
                            tifffile.imwrite(cur_out, merged, shape=merged.shape, imagej=True, photometric='MINISBLACK', metadata={'axes': 'ZYX'})
                        
                    # Store max projection accross the first axis (Z) for current channel 
                    if (self.write_max_projection | self.write_max_projection_onefile):
                        mp = np.max(merged, axis=0)
                        
                        if (self.write_max_projection):
                            log.info(f"{plate} well: {well} channel: {str(channel)} field: {field}, max projection")
                        
                            # Final output filename for channel
                            cur_out = f"{out_dir_final}/{field}_{plate}_{well}_ch{str(channel)}.tiff"
                            
                            if self.uint32:
                                tifffile.imwrite(cur_out, mp, shape=mp.shape, metadata={'axes': 'YX'})        
                            else:
                                tifffile.imwrite(cur_out, mp, shape=mp.shape, imagej=True, metadata={'axes': 'YX'})        
                        else:
                            max_projection.append(mp)

                #-------------------------------------------------
                # Write the max projection one file, CYX
                if (self.write_max_projection_onefile):
                    log.info(f"{plate} well: {well} field: {field} max projection")
                    cur_out = f"{out_dir_final}/{plate}_{well}_{field}_max_projection.tiff"
                    
                    max_merged = np.array(max_projection)
                    
                    if self.uint32:
                        tifffile.imwrite(cur_out, max_merged, shape=max_merged.shape, metadata={'axes': 'CYX'})
                    else:
                        tifffile.imwrite(cur_out, max_merged, shape=max_merged.shape, imagej=True, metadata={'axes': 'CYX'})


# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Stage plate/row/col/field.ome.tiff files into a format compatible for cellprofiler, applying previously calculated registrations, flatifleds and masks")
    parser.add_argument('-w', '--well', help='Well ID to merge', required=True)
    
    # Common argulents
    parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    parser.add_argument('-o','--output', help='Output prefix', required=True)
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+', required=True)
    parser.add_argument('-m','--plate_merge', help='Plates to combine as multiple cycles of the same plate. Output will be seqeuntially added as additional channels in the order provided', nargs='+', default=None)
    parser.add_argument('--registration_dir', help="Path to registration root storing <plate>/<row>/<col>/<field>.pickle", default=None)
    parser.add_argument('--mask_dir', help="Path to mask root storing masks <plate>/<row>/<col>/<field><mask_pattern>.tiff", default=None)
    parser.add_argument('--mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py nucleus. <pattern>.tiff", default="*_nucl_mask_*_cp_masks.tiff")
    parser.add_argument('--mask_channels', help='Channels to mask, output will get 2 extra channels per mask_channel, in the order provided. This is applied after registering and flatfield correction, so use the final channel ids', nargs='+', default=None)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. <plane #1> | [<plane #1> <plane #2> <plane #n>]', nargs='+', default=None)
    parser.add_argument('--flatfields', help='Basicpy model dir for a channel. If merging channels ids are assigned seqeuntially for extra channels <plate>_ch<channel>=/path/to/model ', nargs='+', default=None)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint after applying basicpy model", action='store_true', default=False)

    # Specific arguments
    parser.add_argument('--scaling_factors', help='Divide each channel by this constant factor. Apply this if images are using a small portion of the 16/32 bit space <plate>_ch<channel>=<factor>', nargs='+', default=None)
    parser.add_argument('--no_zstack', help="Do not output one tiff file per channel containing all z-stacks", action='store_true', default=False)
    parser.add_argument('--max_project', help="Output max projection as per channel YX tiffs", action='store_true', default=False)
    parser.add_argument('--max_project_onefile', help="Output max projection one CYX tiff. Specify --no_zstack to skip per channel tiffs.", action='store_true', default=False)
    
    args = parser.parse_args()
  
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    
    runner=MergeAndAlign(args)
    #log.info(f"Well:\t{args.well}")
    #runner.printParams()
    log.info("-----------------------------------------------------------")
    
    runner.run(args.well)








