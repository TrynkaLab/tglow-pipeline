#!/usr/bin/env python

import os
import argparse
import numpy as np

import logging
from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader
from tglow.utils.tglow_utils import float_to_16bit_unint, float_to_32bit_unint, rescale_stack
from aicsimageio.writers import OmeTiffWriter

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def read_and_split_file(filepath, param):
    """
    Reads a file containing a single line of text separated by spaces,
    and returns a list of the items split on the space character.

    Parameters:
        filepath (str): Path to the file.

    Returns:
        list: A list of strings split by space.
    """
    if filepath is None:
        log.warning(f"{param} filepath is none, returning none")
        return None

    if not os.path.isfile(filepath):
        log.warning(f"{param} file not found, returning None. File: {filepath}")
        return None

    with open(filepath, 'r') as file:
        line = file.readline()
        return line.strip().split()


def parse_channel_values(values, plate):
    """Parse a list of '<plate>_ch<channel>=<value>' entries into {channel: value} for one plate."""
    if values is None:
        return None

    prefix = f"{plate}_ch"
    result = {}
    for val in values:
        key, value = val.split("=")
        if key.startswith(prefix):
            channel = int(key[len(prefix):])
            result[channel] = float(value)

    return result


def write_ometiff(stack, plate, row, col, field, output, channel_names, physical_pixel_sizes):
    outdir = f"{output}/{plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}"

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        log.info(f"Folder created: {outdir}")

    OmeTiffWriter.save(stack,
                        f"{outdir}/{field}.ome.tiff",
                        dim_order="CZYX",
                        channel_names=channel_names,
                        physical_pixel_sizes=physical_pixel_sizes)


# Main loop
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Apply per-channel scaling factors to already-finalized unscaled OME-TIFFs. Does not redo flatfield/registration/decon/max-projection.")

    parser.add_argument('-w', '--well', help='Well ID to process', required=True)
    parser.add_argument('-i', '--input', help='Base dir with unscaled OME-TIFFs (plate/row/col/field.ome.tiff)', required=True)
    parser.add_argument('-o', '--output', help='Output base dir for scaled OME-TIFFs', required=True)
    parser.add_argument('-p', '--plate', help='Plate to process', required=True)
    parser.add_argument('--fields', help='Fields to use, defaults to all fields found for the well', nargs='+', default=None)
    parser.add_argument('--pattern', help='Glob pattern used to discover fields', default="*.ome.tiff")
    parser.add_argument('--scaling_factors', help='File with "<plate>_ch<channel>=<factor>" pairs, space separated on one line', required=True)
    parser.add_argument('--scaling_slope', help='File with "<plate>_ch<channel>=<slope>" pairs, controls sigmoid slope. Ignored if --scaling_bias is not set', default=None)
    parser.add_argument('--scaling_bias', help='File with "<plate>_ch<channel>=<bias>" pairs, controls sigmoid bias. Scaling is applied unweighted if this is not set', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint", action='store_true', default=False)
    args = parser.parse_args()

    log.info(f"Input:\t\t{args.input}")
    log.info(f"Output:\t\t{args.output}")
    log.info(f"Plate:\t\t{args.plate}")
    log.info(f"Well:\t\t{args.well}")
    log.info(f"Fields:\t\t{args.fields}")

    factors = parse_channel_values(read_and_split_file(args.scaling_factors, "--scaling_factors"), args.plate)
    slopes = parse_channel_values(read_and_split_file(args.scaling_slope, "--scaling_slope"), args.plate)
    biases = parse_channel_values(read_and_split_file(args.scaling_bias, "--scaling_bias"), args.plate)

    if not factors:
        raise RuntimeError(f"No scaling factors found for plate {args.plate} in {args.scaling_factors}")

    if (biases is None) != (slopes is None):
        raise RuntimeError("Must supply both --scaling_slope and --scaling_bias, or neither")

    reader = AICSImageReader(args.input, plates_filter=[args.plate], fields_filter=args.fields, pattern=args.pattern)

    row, col = ImageQuery.well_id_to_index(args.well)

    for field in reader.fields[args.plate][str(row)][str(col)]:
        log.info(f"Rescaling well {args.well} field {field}")

        iq = ImageQuery(args.plate, row, col, field)
        meta = reader.get_img(iq)
        stack = reader.read_image(iq)

        stack = rescale_stack(stack, factors, slopes, biases)

        if args.uint32:
            stack = float_to_32bit_unint(stack)
        else:
            stack = float_to_16bit_unint(stack)

        write_ometiff(stack, args.plate, row, col, field, args.output,
                       channel_names=meta.channel_names,
                       physical_pixel_sizes=meta.physical_pixel_sizes)

    log.info("Done")
