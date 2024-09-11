#!/usr/bin/env python 

# Copied from Martin Prete's scripts
# pip intstall tifffile imagecodecs ome-model
####################################################################
# Script invocation:
# python actions/companion.py \
#       --input_file="/path/to/stitched/Index.xml" \
#       --output_path="/path/to/place/companion.ome/inside/" 
####################################################################
#
# - The script will crete a companion.ome file from stitched
#   index file after using acapella software.
# - OMERO prefers a companion.ome file to import as a Plate
#   instead of each individual well as an image inside a dataset.
#
####################################################################

import os
import re
import argparse
from pathlib import Path
from tifffile import tiffcomment
from xml.etree import ElementTree as ET
from ome_model.experimental import Plate, Image, create_companion


def main(input_file = "index.xml"):
    images_folder = os.path.dirname(input_file)
    metadata = ET.parse(input_file).getroot()
    default_namespace = re.findall(r'^{(.*)}',metadata.tag)[0]
    NS = {
        "PE": default_namespace,
        "OME": "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    }
    print(f"[*] Reading {input_file}")
    metadata = ET.parse(input_file)

    # collect plate information
    print(f"[*] Reading plate")
    plates = metadata.findall("./PE:Plates/PE:Plate",NS)
    assert len(plates)==1, "Expected only one plate"
    plate = plates[0]
    plate_name = plate.find("./PE:Name",NS).text
    plate_rows = int(plate.find("./PE:PlateRows",NS).text)
    plate_columns = int(plate.find("./PE:PlateColumns",NS).text)
    # gather all images
    images = {}
    for row in range(plate_rows):
        images[str(row+1)] = {}
        for column in range(plate_columns):
            images[str(row+1)][str(column+1)] = {}
    for image in metadata.findall("./PE:Images/PE:Image",NS):
        row = image.find("./PE:Row",NS).text
        column = image.find("./PE:Col",NS).text
        images[row][column] = {
            "filename": f"{images_folder}/{image.find('./PE:URL',NS).text}"
        }
    
    print(f"[*] Building OME Plate")
    # build companion file
    omePlate = Plate(plate_name, plate_rows, plate_columns)
    print(f"[*] Reading image metadata")
    for row in range(plate_rows):
        for column in range(plate_columns):
            well = omePlate.add_well(row, column)
            i = images[str(row+1)][str(column+1)]
            # skip this well if not scanned, meaning it doesn't have matching file ( Images>Image>URL )
            if i.get("filename") is None or not Path(i["filename"]).is_file():
                continue
            well_metadata = ET.fromstring(tiffcomment(i["filename"]))
            pixel_data = well_metadata.find("./OME:Image/OME:Pixels", NS).attrib
            
            image_psx = float(pixel_data['PhysicalSizeX']) if 'PhysicalSizeX' in pixel_data else None
            image_psy = float(pixel_data['PhysicalSizeY']) if 'PhysicalSizeY' in pixel_data else None
            image_psz = float(pixel_data['PhysicalSizeZ']) if 'PhysicalSizeZ' in pixel_data else None
            
            image_name = Path(i['filename']).name
            #print(f"- {image_name}", end="")
            
            image = Image(image_name,
                          sizeX=int(pixel_data['SizeX']),
                          sizeY=int(pixel_data['SizeY']),
                          sizeZ=int(pixel_data['SizeZ']),
                          sizeC=int(pixel_data['SizeC']),
                          sizeT=int(pixel_data['SizeT']),
                          physSizeX=image_psx,
                          physSizeY=image_psy,
                          physSizeZ=image_psz,
                          order=pixel_data['DimensionOrder'],
                          type=pixel_data['Type'])
            
            channels_data = well_metadata.findall("./**/OME:Channel", NS)
            #print(f" - Channels:", end="")
            for channel in channels_data:
                if 'Name' in channel.attrib:
                    image.add_channel(channel.attrib["Name"], channel.attrib["Color"])
                    #print(channel.attrib["Name"], end="|")

            tiffs_data = well_metadata.findall("./**/OME:TiffData", NS)
            for t in tiffs_data:
                tiff_filename = t[0].attrib['FileName'] # <UUID FileName="...">
                image.add_tiff(tiff_filename,
                               c=t.attrib['FirstC'],
                               t=t.attrib['FirstT'],
                               z=t.attrib['FirstZ'],
                               ifd=t.attrib['IFD'],
                               planeCount=t.attrib['PlaneCount'])
            well.add_wellsample(0, image)
            #print()

    companion_filename = f"{images_folder}/{plate_name}.companion.ome"
    print(f"[*] Writing {companion_filename}")
    create_companion(plates=[omePlate], out=companion_filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build OME-XML (companion.ome) from Perkin Elmer index.xml output of stitching')
    parser.add_argument('--input_file', type=str, default="index.xml", required=False, help='Full path to the XML file containing plate information')
    args = parser.parse_args()
    main(args.input_file)

