#!/usr/bin/env python

import tifffile
from PIL import Image
import numpy as np
from scipy.ndimage import center_of_mass
import os
import pandas as pd
import glob
import yaml  # For reading the config file
import h5py
import argparse

import logging
from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def open_tif(path):
    """Open and read an OME-TIFF file, returning the image data as a NumPy array."""
    with tifffile.TiffFile(path) as tif:
        data = tif.asarray()
        # ome_metadata = tif.ome_metadata  # Not used here
    return data

def open_segmentation_mask(path):
    """Open and read a segmentation mask OME-TIFF file."""
    with tifffile.TiffFile(path) as tif:
        mask = tif.asarray()
    return mask



# @author: Francesco Cisterno (wrote originially)
# @author: Olivier Bakker (updated to work in pipeline)
class CellSampler:
    
    def __init__(self, output_root, image_dir, cell_mask_dir, ref_channel, qry_channels, nucl_mask_dir=None, cell_mask_pattern="*_cell_mask_*_cp_masks.tiff", nucl_mask_pattern="*_nucl_mask_*_cp_masks.tiff", image_size=None):
    
        self.output_root=output_root
          
        # Image size is an array of ZYX
        self.image_size=image_size
        self.nuclear_registration = nucl_mask_dir is not None
    
        # Main processing loop      
        self.image_reader = AICSImageReader(image_dir,  pattern="*.ome.tiff")
        self.cell_mask_reader = AICSImageReader(cell_mask_dir,  pattern=cell_mask_pattern)
        
        if nucl_mask_dir is not None:
            self.nucl_mask_reader = AICSImageReader(nucl_mask_dir,  pattern=nucl_mask_pattern)
        else:
            self.nucl_mask_reader = None
        
        if ref_channel is not None:
            self.ref_nucl_channel=int(ref_channel)
            self.qry_nucl_channels=[int(x) for x in qry_channels]
        else:
            self.ref_nucl_channel=None
            self.qry_nucl_channels=[]
        

    def save_cell_crops(self, img, cell_mask, nuclei_mask, iq):
        """
        For each cell in the mask, extract the cell image, compute statistics, and save the result.
        Raw images are saved - so they will not be in the range [0, 255]; moreover, all the channels are saved.
        """
        
        # Create output dir
        os.makedirs(f'{args.output}/{iq.plate}/{iq.get_row_letter()}/{iq.col}/', exist_ok=True)
        
        # Create h5 (one per well) to store cell features
        h5f_path = f'{self.output_root}/{iq.plate}/{iq.get_row_letter()}/{iq.col}/{iq.field}.h5'
        h5f = h5py.File(h5f_path, 'w')
        out_df = pd.DataFrame()
        count = 0
        
        log.info(f"Creating crops over {np.max(cell_mask)} objects")
        for idx in list(range(1, np.max(cell_mask))):  # Skip background (assumed 0)

            # Compute bounding box dimensions (zyx)
            msk_region = np.array(np.where(cell_mask == idx))
            
            # Deals with missing indexes. Not sure why this happens, it seems rare
            if msk_region.size == 0:
                log.warning(f"Masked region is array of size 0, for {idx}. Skipping idx")
                continue
            
            depth = msk_region[1].max() - msk_region[1].min()
            width = msk_region[2].max() - msk_region[2].min()
            height = msk_region[3].max() - msk_region[3].min()
            
            # Deals with single pixel objects by removing them
            if (width == 0) or (height == 0):
                log.warning(f"Masked region has 0 height or width {msk_region.shape} for {idx}. Skipping idx")
                continue
            
            centroid = center_of_mass(cell_mask == idx)
            centroid = [i.round().astype(np.uint32) for i in centroid]
            diameter = max(width, height, depth)

            # Skip cells too close to the border in yx
            if ((centroid[2] - (diameter // 2)) < 0) or ((centroid[2] + (diameter // 2)) > self.image_size[2]) or \
            ((centroid[3] - (diameter // 2)) < 0) or ((centroid[3] + (diameter // 2)) > self.image_size[3]):
                continue
        
            # Crop cell image and masks
            y0 = int(centroid[2] - (diameter // 2))
            y1 = int(centroid[2] + (diameter // 2))
            x0 = int(centroid[3] - (diameter // 2))
            x1 = int(centroid[3] + (diameter // 2))

            # Subset czyx
            cell_img = img[:,:,y0:y1, x0:x1]

            if (cell_img.shape[2] == 0) or (cell_img.shape[3] == 0):
                continue

            cell_mask_bin = ((cell_mask[:,:,y0:y1, x0:x1]) == idx).astype(np.uint8)
            
            #log.info(f"Centroid {centroid}, diameter {diameter}, image size {cell_img.shape}")

            #-------------------------------------------------------------
            # Calculate the correlation between the registration channels
            # Prepare row for output DataFrame
            df_row = {
                'plate': iq.plate, 'row': iq.get_row_letter(), 'col': iq.col, 'field': iq.field, 'well': iq.get_well_id(), 'cell_index': str(idx).zfill(4), 'W': width, 'H': height, 'centroid_z': centroid[1], 'centroid_y': centroid[2], 'centroid_x': centroid[3],  'diameter': diameter
            }

            if nuclei_mask is not None:
                nucl_mask_bin = ((nuclei_mask[:,:,y0:y1, x0:x1]) == idx).astype(np.uint8)

            if self.ref_nucl_channel is not None:
                if self.nuclear_registration:
                    # Compute correlation between two channels in the nucleus region
                    for qry_channel in self.qry_nucl_channels:
                        try:
                            corr = np.corrcoef(cell_img[self.ref_nucl_channel, nucl_mask_bin[0] == cell_mask_bin[0]], cell_img[qry_channel, nucl_mask_bin[0] == cell_mask_bin[0]])[1, 0]
                        except Exception:
                            corr = np.nan
                            
                        df_row[f"corr_{self.ref_nucl_channel}_{qry_channel}"] = corr
                else:
                    # If there is no nucleus, calculate it in the cell region
                    for qry_channel in self.qry_nucl_channels:
                        try:
                            corr = np.corrcoef(cell_img[self.ref_nucl_channel,cell_mask_bin[0]], cell_img[qry_channel, cell_mask_bin[0]])[1, 0]
                        except Exception:
                            corr = np.nan         
                        
                        df_row[f"corr_{self.ref_nucl_channel}_{qry_channel}"] = corr
                        
                    
            #-------------------------------------------------------------
            # Get max value per channel
            cell_channel_max = np.max(cell_img, axis=(1,2,3))
                    
            for ch_idx, x in enumerate(cell_channel_max):
                df_row[f"max{ch_idx}"] = x
            
            #-------------------------------------------------------------            
            # Add cell and nuclei segmentation masks as last channel
            cell_img = np.concatenate((cell_img, cell_mask_bin), 0)
            
            if nuclei_mask is not None:
                cell_img = np.concatenate((cell_img, nucl_mask_bin), 0)

            #log.info(f"Final image {idx} of shape {cell_img.shape} with corr {corr}")
            
            # Save cell image as dataset in the group /row/col/image
            h5f.create_dataset(str(idx).zfill(4), data=cell_img, compression='gzip')

            # Append row to DataFrame
            out_df = pd.concat((out_df, pd.DataFrame(df_row, index=[count])))
            count += 1
            
        out_df.to_csv(f'{self.output_root}/{iq.plate}/{iq.get_row_letter()}/{iq.col}/{iq.field}.csv')
        
        # Add channel names
        img_meta = self.image_reader.get_img(iq)
        
        channel_names = img_meta.channel_names #list(range(0,img.shape[0]))
        channel_names.append("cell_mask")
        if nuclei_mask is not None:
            channel_names.append("nucl_mask")
        channel_names = [str(x) for x in channel_names]
        
        h5f.create_dataset("channel_names", data=channel_names)
        
        h5f.flush()
        h5f.close()
        
        # Experimented with this
        #out_df.to_hdf(h5f_path, key='/cell_index', mode='a')
        
    def sample_imageset(self, iq, max_cells_per_field=1000):
        
        # Read image CZYX
        img = self.image_reader.read_stack(iq) 

        # If the image size is not set, set it here
        if self.image_size is None:
            self.image_size = img.shape

        # Load masks
        mask = self.cell_mask_reader.read_stack(iq)

        # Get cells number
        cells_number = np.max(mask)
        
        log.info(f"Read cellmask. It has {cells_number} objects in it")
        
        if cells_number > max_cells_per_field:
            return None

        # Nuclei mask
        if self.nucl_mask_reader is not None:
            nucl_mask = self.nucl_mask_reader.read_stack(iq)
            # Only keep nuclei inside cells, and re-assign to the cell label
            nucl_mask = ((nucl_mask > 0) & (mask > 0)) * mask
        else:
            nucl_mask = None
            
        log.info(f"Read images of shape {img.shape}, cell_masks {mask.shape}, nucl_masks {nucl_mask.shape}")
        
        self.save_cell_crops(img, mask, nucl_mask, iq)

    
if __name__ == "__main__":

    # Create the parser
    parser = argparse.ArgumentParser(description="Produce cellcrops from /plate/row/col/field.ome.tiff images and /plate/row/col/field.ome.tiff masks. Must be CZYX")

    # Add arguments
    parser.add_argument('-i','--input', help='Base dir to input images CYZX ome tiffs', required=True)
    parser.add_argument("--plate", help="Experiment plate", required=True)
    parser.add_argument("--well", help="Image well", required=True)
    parser.add_argument('-o','--output', help='Output prefix', required=True)
    parser.add_argument('--cell_mask_dir', help="Path to mask root storing masks <plate>/<row>/<col>/<field><mask_pattern>.tiff", required=True)
    parser.add_argument('--cell_mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py nucleus. <pattern>.tiff", default="*_cell_mask_*_cp_masks.tiff")
    parser.add_argument('--nucl_mask_dir', help="Path to mask root storing masks <plate>/<row>/<col>/<field><mask_pattern>.tiff", default=None)
    parser.add_argument('--nucl_mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py nucleus. <pattern>.tiff", default="*_nucl_mask_*_cp_masks.tiff")
    parser.add_argument('--ref_channel', help="The registration channel in reference plate", required=False, default=None)
    parser.add_argument('--qry_channels', help="The registration channel in query plates", required=False, nargs="+", default=[])
    parser.add_argument('--max_per_field', help="Max number of cells per field to consider", default=5000)

    # Parse the arguments
    args = parser.parse_args()

    # Get plate, row and column
    iq = ImageQuery.from_plate_well(args.plate, args.well)

    sampler = CellSampler(args.output,args.input, args.cell_mask_dir, args.ref_channel, args.qry_channels, args.nucl_mask_dir,  args.cell_mask_pattern, args.nucl_mask_pattern)

    for imagename in glob.glob(os.path.join(args.input, iq.plate, iq.get_row_letter(), iq.col, '*.ome.tiff')):
        # Get image id
        field = imagename.split('/')[-1].split('.ome.tiff')[0].split('/')[-1]
        iq.field = field
        log.info(f"iq: {iq.to_string()}")
        sampler.sample_imageset(iq, int(args.max_per_field))
        
    print('Completed without errors', flush=True)
