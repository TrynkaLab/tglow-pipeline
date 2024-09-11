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
import pystackreg
import pickle
import copy
import warnings
import csv
import xml.etree.ElementTree as ET

#-------------------------------------------------------------------
# Static functions
#-------------------------------------------------------------------

# Plot before and after of an imageset for a given index in the arrays
def plot_imgs(before, after, filename):
    
    #print(f"[DEBUG] b:{before.dtype} a: {after.dtype}")
    
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow((before * 255).astype(np.uint8), cmap='gray', vmin=0, vmax=1)
    axes[0].set_title("Before registration")
    im = axes[1].imshow((after * 255).astype(np.uint8), cmap='gray', vmin=0, vmax=1)
    axes[1].set_title("After registration")
    
    fig.tight_layout()
    fig.savefig(filename, dpi=600)
    plt.close()

def composite_images(imgs, equalize=False, aggregator=np.mean):

    if equalize:
        imgs = [exposure.equalize_hist(img) for img in imgs]

    imgs = [img / img.max() for img in imgs]

    if len(imgs) < 3:
        imgs = [np.zeros(shape=imgs[0].shape)] * (3-len(imgs)) + imgs

    imgs = np.dstack(imgs)
    return imgs

def model_dict_to_str(dict):
    
    if not dict == None:
        str="{"
        for key in dict.keys():
            str=str+f"'{key}': {type(dict[key])},"
        str=str+"}"
        return str
    else:
        return None
    
# Build dict linking letters to number
def build_well_index(upper=True):
    if (upper):
        chars = string.ascii_uppercase
    else:
        chars = string.ascii_lowercase
    nums = list(range(1, len(chars)))
    idx = dict(zip(chars,nums))
    return(idx)

# Function that converts wells 
def make_prefix(well, idx, prepend_zero=True):
    match = re.search("([a-zA-Z])(\\d+)", well)
    if match != None:
        row = match[1]
        col = int(match[2])
        row_idx = idx[row]
        if (prepend_zero):
            prefix = "r" + str(row_idx).zfill(2) + "c" + str(col).zfill(2)
        else:
            prefix = "r" + str(row_idx) + "c" + str(row)
        return prefix
    else:
        return None

# Encode numpy formats as JSON
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

# Get channel info easliy queriable in python as json format
# adapted from Matin Prete's ParseXML.py
def get_channel_channel_info(index_xml) -> dict:
    tree = ET.parse(index_xml)
    root = tree.getroot()

    default_namespace = re.findall(r'^{(.*)}',root.tag)[0]
    NS = {
        "PE": default_namespace,
        "OME": "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    }

    entries=root.findall("PE:Maps/PE:Map/PE:Entry", NS)

    # Adapted & extended from Martins script
    channel_info={}
    for e in entries:
        cn = e.findall('./PE:ChannelName',NS)
        
        if cn:
            channel = {
                "id": e.attrib['ChannelID'],
                "name": cn[0].text,
                "image_type": e.find('./PE:ImageType',NS).text,
                "acquisition_type": e.find('./PE:AcquisitionType',NS).text,
                "illumination_type": e.find('./PE:IlluminationType',NS).text,
                "channel_type": e.find('./PE:ChannelType',NS).text,
                "binning_x": int(e.find('./PE:BinningX',NS).text),
                "binning_y": int(e.find('./PE:BinningY',NS).text),
                "emmisson": int(e.find('./PE:MainEmissionWavelength',NS).text),
                "excitation": int(e.find('./PE:MainExcitationWavelength',NS).text),
                "max_intensity": int(e.find('./PE:MaxIntensity',NS).text),
                "numerical_apeture": float(e.find('./PE:ObjectiveNA',NS).text),
                "objective_magnification": float(e.find('./PE:ObjectiveMagnification',NS).text),
                "exposure_time":{"unit": e.find('./PE:ExposureTime',NS).attrib['Unit'],
                                "value": float(e.find('./PE:ExposureTime',NS).text)},
                "size": (
                    int(e.find('./PE:ImageSizeX',NS).text),
                    int(e.find('./PE:ImageSizeY',NS).text)
                )}
            
            x = e.find('./PE:ImageResolutionX',NS)
            y = e.find('./PE:ImageResolutionY',NS)
            channel["image_resolution"] = {
                "x": {"unit":x.attrib["Unit"], "value":float(x.text),},
                "y": {"unit":y.attrib["Unit"], "value":float(y.text),},
            }
            # more info inside :
            ff_selector = f"./PE:Maps/PE:Map/PE:Entry[@ChannelID='{channel['id']}']/PE:FlatfieldProfile"
            channel["flatfield"] = root.find(ff_selector,NS).text
            channel_info[channel['id']] = channel
            
    return channel_info

# Convert a numpy 32bit float to a 16 bit int, clip values to zero and 65535
def float_to_16bit_unint(matrix) -> np.array:

    # Convert to 16 bit to keep consistency with output
    # Round to nearest int            
    matrix = np.rint(matrix)
    
    # Set negatives to zero
    matrix[matrix < 0] = 0
    
    # Clip values at the 16 bit max
    matrix[matrix > np.iinfo(np.uint16).max] = np.iinfo(np.uint16).max
    
    matrix = matrix.astype(np.uint16)
    
    return matrix

# Convert a numpy 32bit float to a 16 bit int, clip values to 
def float_to_32bit_unint(matrix) -> np.array:

    # Round to nearest int            
    matrix = np.rint(matrix)
    
    # Set negatives to zero
    matrix[matrix < 0] = 0
    
    # Clip values at the 32 bit max
    matrix[matrix > np.iinfo(np.uint32).max] = np.iinfo(np.uint32).max
    
    matrix = matrix.astype(np.uint32)
    
    return matrix


# Read an array of tiff files as a 3d numpy array
# Optionally transform them using a precomputed alignment matrix
def read_tiffstack_as_numpy(cur_files, alignment=None):
    i=0
    arrays=[]    
    # Loop over the raw tiffs and read them as 2d numpy arrays
    while i < len(cur_files):      
        cur_img = tifffile.imread(cur_files[i])
        if not alignment is None:
            #sr = StackReg(transform)
            #cur_img= sr.transform(cur_img, tmat=alignment)
            tform = transform.AffineTransform(matrix=alignment)

            # transform image using the saved transformation matrix
            cur_img = transform.warp(cur_img, tform, preserve_range=True)
                
        arrays.append(cur_img)
        i=i+1
    
    # Merge the array of 2d numpy arrays into a 3d array
    merged=np.array(arrays)
    
    return(merged)


#-------------------------------------------------------------------
# Main runner
class MergeAndAlign:
    
    COLLUMN_INDEX=build_well_index()
    FIELD_PREFIX="f"
    PLANE_PREFIX="p"
    CHANN_PREFIX="ch"
    
    # Constructor, and parse arguments
    def __init__(self, args):
          
        self.input=args.input
        self.output=args.output
        self.write_zstack=not args.no_zstack
        self.write_max_projection=args.max_project
        self.uint32=args.uint32
        
        # Set subfolders / plates to process
        if args.plate == None:
            self.plates=[ name for name in os.listdir(input) if os.path.isdir(os.path.join(input, name)) ]
        else:
            self.plates=args.plate

        
        # Channels for alignment and plates to merge
        self.plates_merge=args.plate_merge   
        if not self.plates_merge is None:
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
            self.ref_channel=None
            self.qry_channel=None
            self.ref_channel_eval=None
            self.qry_channel_eval=None

        # Fields, planes, channels
        self.fields=args.fields
        self.planes=args.planes
        self.channels=[int(x) for x in args.channels]
        
        if not args.channels_merge is None:
            self.channels_merge=[int(x) for x in args.channels_merge]
        else:
            self.channels_merge=None

        # Registration options
        self.plot=True
        self.save_reg=False
        self.transform=StackReg.TRANSLATION

        # Build dict with pretrained BaSiCpy models
        if args.basicpy_model == None:
            self.flatfields=None
        else:
            self.flatfields={}
            for val in args.basicpy_model:
                keypair = val.split("=")
                print(f"[INFO] {keypair}")
                #basicpy_models[keypair[0]] = load_model(keypair[1])
                self.flatfields[keypair[0]] = BaSiC.load_model(keypair[1])
                self.flatfields[keypair[0]].baseline="" # Hack arround missing baseline    
                
        self.eq_merge=args.eq_merge   
        self.inv_merge=args.inv_merge   
        self.eval_merge=args.eval_merge
        
        if self.eval_merge:
            if (args.ref_channel_eval is None) | (args.qry_channel_eval is None):
                raise RuntimeError("Need to supply reference and query evaluation channels for registration evaluation")
            

    # Main loop
    def run(self, well):
        # Construct prefix for well
        self.prefix=make_prefix(well, MergeAndAlign.COLLUMN_INDEX)

        # Index with image files 
        channels=self.channels.copy()
        index = self.build_plate_index(well, self.plates, self.channels)
        
        # If multiple channels should be merged, run the registration routine
        final_alignment=None
        channel_index_new=None
        
        if (self.plates_merge != None):
            
            if (len(index) !=1):
                raise RuntimeError("When merging channels --plate can only have one value as it is used as the reference plate")
            
            # Build the index for the plates to register and merge
            index_merge = self.build_plate_index(well, self.plates_merge, self.channels_merge)
        
            ref_plate = self.plates[0]
            alignment_matrices = self.calculate_registration(index, index_merge, well)
            
            # Merge indexes to treat the extra plates as channels of the same plate
            # This ensures the downstream code still works, and produces tiff stacks as if the images
            # came from the same plate
            index_new = copy.deepcopy(index)
            channels_new = self.channels.copy()
            channel_start = np.max(channels)
            cycle=['cycle_1' for ch in channels]
            channel_index_new={}
                    
            for plate_merge in self.plates_merge:
                i=0;
                channel_index_new[plate_merge] = {}
                for field in self.fields:
                    field=MergeAndAlign.FIELD_PREFIX + str(field).zfill(2)
                    newch = channel_start + 1

                    for channel_merge in self.channels_merge:
                        if i==0:
                            channel_index_new[plate_merge][channel_merge]=newch
                            channels_new.append(newch)
                            cycle.append(plate_merge)
                        index_new[ref_plate][well][field][newch] = index_merge[plate_merge][well][field][channel_merge]
                        newch+=1
                        
                    i += 1
                channel_start = np.max(channels_new)
                
            # Re-code the alignment to the new channel indices per channel thing
            # Channels in the original plate get set to None
            final_alignment = {}
            final_alignment[ref_plate] = {}
            final_alignment[ref_plate][well] = {}
            final_alignment['transform'] = self.transform
            
            for field in self.fields:
                field=MergeAndAlign.FIELD_PREFIX + str(field).zfill(2)
                final_alignment[ref_plate][well][field] = {}
                
                i=0
                while i < len(channels_new):
                    channel = channels_new[i]
                    cur_cycle = cycle[i]
                    if len(final_alignment[ref_plate][well][field]) == 0:
                        final_alignment[ref_plate][well][field][channel] = {}
                        
                    if cur_cycle == 'cycle_1':
                        final_alignment[ref_plate][well][field][channel] = None
                    else:
                        final_alignment[ref_plate][well][field][channel] = alignment_matrices[cur_cycle][field]
                    i+=1

            # Update CHANNELS, index
            channels=channels_new
            index=index_new
            
        print("[INFO] Writing channel info")
        self.write_updated_channel_indices(channel_index_new)

        if (not (self.write_max_projection | self.write_zstack)):
            print("[INFO] No output specified, exiting")
            return


        print("[INFO] Staring merging")
        
        if final_alignment != None:
            print("[INFO] Applying precalculated transformations while reading")

        # Loop over possible plates
        for plate in self.plates:
            out_dir_final= f"{self.output}/{plate}/{well}/"
            
            if not os.path.exists(out_dir_final):
                os.makedirs(out_dir_final)
                print(f"[INFO] Folder created: {out_dir_final}")
                    
            for field in self.fields:
                
                field=MergeAndAlign.FIELD_PREFIX + str(field).zfill(2)

                # Array to store the fields for the max projection
                max_projection = []
                
                # Loop over channels
                for channel in channels:
                
                    print(f"[INFO] {plate} well: {well} channel: {str(channel)} field: {field}")
                    
                    # Final output filename for channel
                    cur_out = f"{out_dir_final}/{plate}_{well}_{field}_ch{str(channel)}.tiff"      
                    
                    if (final_alignment != None):
                        merged = read_tiffstack_as_numpy(index[plate][well][field][channel],
                                                         final_alignment[plate][well][field][channel])
                        merged = float_to_16bit_unint(merged)
                    else:
                        merged = read_tiffstack_as_numpy(index[plate][well][field][channel])

                    # If flatfields are specified, apply them here
                    if not self.flatfields == None:
                        if str(channel) in self.flatfields.keys():
                            print(f"[INFO] Applying BaSiC model for channel: {str(channel)}")
                            basic_model = self.flatfields[str(channel)]
                            merged = basic_model.transform(merged)
                            
                            # Convert to either 32 or 16 bit unsigned int
                            if self.uint32:
                                merged=float_to_32bit_unint(merged)
                            else:
                                merged=float_to_16bit_unint(merged)
                                
                    # Write out the final tiff file, saves as ZYX
                    if self.write_zstack:
                        tifffile.imwrite(cur_out, merged, shape=merged.shape, imagej=True, photometric='MINISBLACK', metadata={'spacing': 0.149, 'unit': 'um','axes': 'ZYX'})
                        
                    # Store max projection accross the first axis (Z) for current channel 
                    if (self.write_max_projection):
                        mp = np.max(merged, axis=0)
                        max_projection.append(mp)

                # Write the max projection
                if (self.write_max_projection):
                    print(f"[INFO] {plate} well: {well} field: {field} max projection")
                    cur_out = f"{out_dir_final}/{plate}_{well}_{field}_max_projection.tiff"
                    
                    max_merged = np.array(max_projection)
                    tifffile.imwrite(cur_out, max_merged, shape=max_merged.shape, imagej=True, metadata={'spacing': 0.149, 'unit': 'um','axes': 'CYX'})

    # Read PE index file for Phenix or Operetta
    def read_index_file(self, plate):
        
        cur = None
        # Phenix and operetta call their metadata files differently
        if (os.path.isfile(f"{self.input}/{plate}/Index.xml")):
            cur = get_channel_channel_info(f"{self.input}/{plate}/Index.xml")
        elif (os.path.isfile(f"{self.input}/{plate}/Index.idx.xml")):
            cur = get_channel_channel_info(f"{self.input}/{plate}/Index.idx.xml")
        else:
            Warning("[WARN] index xml file not found")
         
        return(cur)


    # Perform registration between planes
    def calculate_registration(self, index, index_merge, well):
            
        alignment_matrices={}

        # Construct the output dir
        plate = self.plates[0]
        out_dir_final=f"{self.output}/{plate}/{well}/"
        reg_dir=f"{out_dir_final}/registration/"

        if not os.path.exists(reg_dir):
            os.makedirs(reg_dir)
            print(f"[INFO] Folder created: {reg_dir}")
        
        if len(index_merge[self.plates_merge[0]]) != len(index[plate]):
            raise RuntimeError("Different number of fields for qry and ref plate")
        
        if self.eval_merge:
            outfile = f"{reg_dir}/{well}_registration_eval_pearson_correlations.tsv"
            file = open(outfile, "w")
            file.write("image\tref_plate\tqry_plate\tref_ch_eval\tqry_ch_eval\tref_ch\tqry_ch\tr_before\tp_before\tr_after\tt_after\n")
            file.close()
            
        print("[INFO] Running alignment between plates on max projections")
        i=0; 
        
        for field in self.fields:
            field = MergeAndAlign.FIELD_PREFIX + str(field).zfill(2)
            
            ref_mp = np.max(read_tiffstack_as_numpy(index[plate][well][field][self.ref_channel]), axis=0)   
            
            
            
            # Loop over plates to merge onto ref
            for plate_merge in self.plates_merge:
                print(f"[INFO] Processing field {field} plate {plate_merge}")
                if i==0:
                    alignment_matrices[plate_merge] = {} 
                    i+=1
                
                qry_mp = np.max(read_tiffstack_as_numpy(index_merge[plate_merge][well][field][self.qry_channel[plate_merge]]), axis=0)
                
                # Run histogram matching between query and reference
                if self.eq_merge:
                    qry_mp = exposure.match_histograms(qry_mp, ref_mp)
                    qry_mp = float_to_16bit_unint(qry_mp)
                    
                # Invert the images priror to registration
                if self.inv_merge:
                    
                    #print(f"[INFO] q:{qry_mp.dtype} r:{ref_mp.dtype}")
                    #qry_mp = np.iinfo(qry_mp.dtype).max - qry_mp
                    #ref_mp_final = np.iinfo(ref_mp.dtype).max - qry_mp
                    ref_mp_final = exposure.rescale_intensity(ref_mp)
                    qry_mp = exposure.rescale_intensity(qry_mp, in_range=(ref_mp.min(), ref_mp.max()))
                    
                    qry_mp = util.invert(qry_mp)
                    ref_mp_final = util.invert(ref_mp_final)
                    
                    thresh = filters.threshold_otsu(ref_mp_final)
                    print(f"[INFO] Post inversion otsu thresholds {thresh}")
                    ref_mp_final[ref_mp_final < thresh] = 0
                    qry_mp[qry_mp < thresh] = 0
                    #ref_mp_final[ref_mp_final > thresh[1] & ref_mp_final < thresh[2]] = 0
                    #qry_mp[qry_mp > thresh[1] & qry_mp < thresh[2]] = 0
                    #print(f"[INFO] q:{qry_mp.dtype} r:{ref_mp.dtype}")

                else:
                    ref_mp_final = ref_mp
               
                # Calculate the offsets
                sr = StackReg(self.transform)
                align_mat = sr.register(ref_mp_final, qry_mp)
                                
                # Save the offsets
                alignment_matrices[plate_merge][field] = align_mat

                # Plot 
                if self.plot:
                    # Apply the offesets
                    qry_mp_reg = sr.transform(qry_mp, tmat=align_mat)
                    before_reg = composite_images([ref_mp_final, qry_mp])
                    after_reg = composite_images([ref_mp_final, qry_mp_reg])
                    plot_imgs(before_reg, after_reg, filename=f"{reg_dir}/{plate_merge}_{well}_{field}_refch{self.ref_channel}_qrych{self.qry_channel[plate_merge]}.png")
                
                if self.eval_merge:
                    print(f"[INFO] Calculation pearson correlation between ref and qry for channels {self.ref_channel_eval} and {self.qry_channel_eval[plate_merge]}")
                    a = np.max(read_tiffstack_as_numpy(index[plate][well][field][self.ref_channel_eval]), axis=0)   
                    b = np.max(read_tiffstack_as_numpy(index_merge[plate_merge][well][field][self.qry_channel_eval[plate_merge]], alignment=align_mat), axis=0)   
                    c = np.max(read_tiffstack_as_numpy(index_merge[plate_merge][well][field][self.qry_channel_eval[plate_merge]]), axis=0)   

                    
                    apcc, apval = measure.pearson_corr_coeff(a, c)
                    print(f"[INFO] Before Pearson R: {apcc:0.3g}, p-val: {apval:0.3g}")
    
                    bpcc, bpval = measure.pearson_corr_coeff(a, b)
                    print(f"[INFO] After Pearson R {bpcc:0.3g}, p-val: {bpval:0.3g}")
                    
                    line=f"{well}_{field}\t{plate}\t{plate_merge}\t{self.ref_channel_eval}\t{self.qry_channel_eval[plate_merge]}\t{self.ref_channel}\t{self.qry_channel[plate_merge]}\t{apcc}\t{apval}\t{bpcc}\t{bpval}\n"
                    file = open(outfile, 'a')
                    file.write(line)
                    file.close()
                    
                    
        if self.save_reg:
            alignment_matrices['transform'] = transform
            pickle.dump(alignment_matrices, open(f"{reg_dir}/registration_{well}.pickle", "wb")) 
            
        return alignment_matrices

    # Write a json file with the channel id's in the merged tiffs
    def write_updated_channel_indices(self, channel_index_new):
        
        channel_info = {}
        for plate in self.plates:
            cur = self.read_index_file(plate)
            
            if (cur != None):
                for key in cur.keys():
                    cur[key]['cycle']=1
                    cur[key]['plate']=plate
                channel_info = {**channel_info, **cur}
            else:
                Warning(f"[WARN] Index file for plate {plate} not found. Skipping writing index")
                return None

        if not channel_index_new is None:
            i=2
            for plate in self.plates_merge:
                cur = self.read_index_file(plate)
                if (cur != None):
                    for key in cur.keys():
                        cur[key]['cycle']=i
                        cur[key]['plate']=plate
                        cur[key]['id_old'] = cur[key]['id']
                        cur[key]['id'] = channel_index_new[plate][int(key)]
                        channel_info[channel_index_new[plate][int(key)]] = cur[key]
                    i+=1
                else:
                    Warning(f"[WARN] Index file for merge plate {plate} not found. Skipping writing index")
                    return None 
                   
        # Convert ids to string for JSON compatability
        channel_info = {str(k):v for k,v in channel_info.items()}           
                
        dir=f"{self.output}/{self.plates[0]}/"
        if not os.path.exists(dir):
            os.makedirs(dir)
            print(f"[INFO] Folder created: {dir}")
        
        file=open(f"{dir}/channel_info.json", "w") 
        json.dump(channel_info, file, indent=4,  cls=NpEncoder)
        file.close()

    # Returns a dict structure for the items in a well and plate 
    # {plate}/{well}/{field}/{channel}/[plane]
    # plane is ordered by plane number
    # TODO: do this based on the PE index xml, should be pretty easy and a cleaner solution
    def build_plate_index(self, well, plates, channels):
        index = {}
        i=0
        for plate in plates:
            cur_path = f"{self.input}/{plate}/{well}/"
            files = os.listdir(cur_path)
            tmp_plate = {}
            
            if i == 0:
                index[plate] = {}
                index[plate][well]={}
                
            for field in self.fields:
                field=MergeAndAlign.FIELD_PREFIX + str(field).zfill(2)
                
                tmp_field={}
                for channel in channels:
                    search_pattern = self.prefix + field + MergeAndAlign.PLANE_PREFIX + "\\d\\d" + "-" + MergeAndAlign.CHANN_PREFIX + str(channel) + "sk1fk1fl1.tiff"
                    # List all the files in that folder
                    cur_files = [file for file in files if re.search(search_pattern, file) != None]
                    cur_planes = [int(re.search(self.prefix + field + MergeAndAlign.PLANE_PREFIX + "(\\d+)" + "-.*", file).group(1)) for file in cur_files]
                    
                    # Find the files that are relevant
                    cur_files = [x for _, x in sorted(zip(cur_planes, cur_files))]
                    cur_files = [cur_path + file for file in cur_files]
                    
                    tmp_field[channel] = cur_files
                    
                    tmp_plate[field] = tmp_field
            
            index[plate][well] = tmp_plate
        return index
    
    def printParams(self):
        
        print(f"Input plates:\t{str(self.plates)}")
        print(f"Input fields:\t{str(self.fields)}")
        print(f"Input planes:\t{str(self.planes)}")   
        print(f"Input channel:\t{str(self.channels)}")     
        print(f"Input:\t\t{self.input}")
        print(f"Output:\t\t{self.output}")    
        print(f"Basicpy models:\t{model_dict_to_str(self.flatfields)}")    
        print(f"Output 32bit:\t{str(self.uint32)}")    
        
        print(f"----------------------------")    
        print(f"Merge plates:\t{str(self.plates_merge)}")     
        print(f"Merge channels:\t{str(self.channels_merge)}")  
        print(f"Ref channel:\t{str(self.ref_channel)}")  
        print(f"Qry channel:\t{str(self.qry_channel)}")
        print(f"Ref channel eval:\t{str(self.ref_channel_eval)}")  
        print(f"Qry channel eval:\t{str(self.qry_channel_eval)}")  
        print(f"Merge hist match:\t{str(self.eq_merge)}")  
        print(f"Merge invert:\t{str(self.inv_merge)}")  
        print(f"Transform:\t{str(self.transform)}")        
        print(f"Save registr:\t{str(self.save_reg)}")     
        print(f"Plot:\t\t{str(self.plot)}")     

# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Merge raw Phenix tiff files organized per well into per channel tiffs, one page per z-stack")
    parser.add_argument('-w', '--well', help='Well ID to merge')
    parser.add_argument('-i','--input', help='Base dir to raw input')
    parser.add_argument('-o','--output', help='Output prefix')
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+')
    parser.add_argument('-r','--ref_channel', help='Reference channel for registration (nucleus)', nargs=1, default=None)
    parser.add_argument('-q','--qry_channel', help='Query channel for registration (nucleus)', nargs=1, default=None)
    parser.add_argument('-m','--plate_merge', help='Plates to align using pystackreg (e.g. multiple cycles of the same plate). Output will be added as additional channels.', nargs='+', default=None)
    parser.add_argument('--no_zstack', help="Do not output one tiff file per channel containing all z-stacks", action='store_true', default=False)
    parser.add_argument('--max_project', help="Output max projection in parallel to per channel tiffs", action='store_true', default=False)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=[1,2,3,4,5,6,7,8,9])
    parser.add_argument('--planes', help='Z planes to use. <plane #1> | [<plane #1> <plane #2> <plane #n>]', nargs='+', default=[1,2,3,4,5])
    parser.add_argument('--channels', help='Channels to use. <channel #1> | [<channel #1> <channel #2> <channel #n>]', nargs='+', default=[1,2,3,4,5])
    parser.add_argument('--channels_merge', help='Channels to for merging from second plate. <channel #1> | [<channel #1> <channel #2> <channel #n>]', nargs='+', default=None)
    parser.add_argument('--basicpy_model', help='Basicpy model dir for a channel. If merging channels ids are assigned seqeuntially for extra channels <channel #1>=/path/to/model | [<channel #1>=/path/to/model <channel #n>=/path/to/model]', nargs='+', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint after applying basicpy model", action='store_true', default=False)
    parser.add_argument('--eq_merge', help="Run histogram matching prior to registration.", action='store_true', default=False)
    parser.add_argument('--inv_merge', help="Invert image prior to registration", action='store_true', default=False)
    parser.add_argument('--eval_merge', help="Calculate correlation between two channels before and after", action='store_true', default=False)
    parser.add_argument('--ref_channel_eval', help='Reference channel for registration evaluation (nucleus)', nargs=1, default=None)
    parser.add_argument('--qry_channel_eval', help='Query channel for registration evaluation (nucleus)', nargs=1, default=None)
    args = parser.parse_args()
  
    print("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        print(f"Folder created: {args.output}")
    

    runner=MergeAndAlign(args)
    print(f"Well:\t{args.well}")
    runner.printParams()
    print("-----------------------------------------------------------")
    
    runner.run(args.well)

    #main(well, input, output, ref_channel, qry_channel, write_zstack, write_max_projection, basicpy_models, uint32)






