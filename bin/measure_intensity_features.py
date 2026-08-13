#!/usr/bin/env python

import os
import argparse
import random
import pandas as pd

import logging
from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# DUMMY IMPLEMENTATION: real object detection / intensity computation is not implemented yet.
# All feature values are populated with this placeholder instead.
DUMMY_VALUE = 1.0

# DUMMY IMPLEMENTATION: per-cell registration correlation isn't computed yet either -
# fabricate one ch1_ch<N>__registration_corr column per non-reference channel (2..N),
# with a spread of values either side of the QC report's default qc_regcor=0.6, so the
# registration QC tab has something meaningful to filter/plot pending the real thing.
DUMMY_REGISTRATION_CORR_RANGE = (0.3, 1.0)

# DUMMY IMPLEMENTATION: no masks are read yet, so this many placeholder objects are
# fabricated per field instead of being counted from a real cell mask.
N_DUMMY_OBJECTS_PER_FIELD = 100

OBJECT_STATS = ["min", "q25", "q50", "mean", "q75", "q99", "q99.9", "q99.999", "max"]
IMAGE_STATS = OBJECT_STATS + ["background", "otsu_raw", "otsu_cellrm", "otsu_artefactrm"]


def channel_columns(n_channels, stats):
    """Build the ch1__<stat>, ch2__<stat>, ... column names for the given number of channels."""
    columns = {}
    for channel in range(1, n_channels + 1):
        for stat in stats:
            columns[f"ch{channel}__{stat}"] = DUMMY_VALUE
    return columns


def registration_corr_columns(n_channels):
    """DUMMY: one ch1_ch<N>__registration_corr value per non-reference channel."""
    return {
        f"ch1_ch{channel}__registration_corr": random.uniform(*DUMMY_REGISTRATION_CORR_RANGE)
        for channel in range(2, n_channels + 1)
    }


def get_n_channels(reader, plate, row, col, fields):
    """Determine the channel count from the first field found for this well."""
    if not fields:
        return 0

    iq = ImageQuery(plate, row, col, fields[0])
    return len(reader.get_img(iq).channel_names)


def main(args):

    well_iq = ImageQuery.from_plate_well(args.plate, args.well)
    row = well_iq.row
    col = well_iq.col
    row_letter = well_iq.get_row_letter()
    well_id = well_iq.get_well_id()

    reader = AICSImageReader(args.input, plates_filter=[args.plate])
    fields = reader.fields.get(args.plate, {}).get(row, {}).get(col, [])

    n_channels = get_n_channels(reader, args.plate, row, col, fields)

    log.info(f"Plate {args.plate} well {args.well}: found {len(fields)} fields, {n_channels} channels")

    object_rows = []
    image_rows = []

    for field in fields:

        # Image-level row - one per field
        image_rows.append({
            "plate": args.plate,
            "row": row_letter,
            "col": col,
            "field": field,
            "well": well_id,
            **channel_columns(n_channels, IMAGE_STATS)
        })

        # Object-level rows - N_DUMMY_OBJECTS_PER_FIELD placeholder objects per field
        for obj_idx in range(1, N_DUMMY_OBJECTS_PER_FIELD + 1):
            object_rows.append({
                "object_id": f"{field}_{str(obj_idx).zfill(4)}",
                "plate": args.plate,
                "row": row_letter,
                "col": col,
                "field": field,
                "well": well_id,
                **channel_columns(n_channels, OBJECT_STATS),
                **registration_corr_columns(n_channels)
            })

    object_df = pd.DataFrame(object_rows)
    image_df = pd.DataFrame(image_rows)

    outdir = os.path.join(args.output, args.plate, row_letter, col)
    os.makedirs(outdir, exist_ok=True)

    object_path = os.path.join(outdir, "object_features.parquet")
    image_path = os.path.join(outdir, "image_features.parquet")

    object_df.to_parquet(object_path, index=False)
    image_df.to_parquet(image_path, index=False)

    log.info(f"Wrote {object_path} ({len(object_df)} rows)")
    log.info(f"Wrote {image_path} ({len(image_df)} rows)")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="DUMMY: produce placeholder per-object and per-image intensity feature tables for a well.")

    parser.add_argument('-i', '--input', help='Base dir to input images plate/row/col/field.ome.tiff', required=True)
    parser.add_argument('-p', '--plate', help='Plate to process', required=True)
    parser.add_argument('-w', '--well', help='Well ID to process', required=True)
    parser.add_argument('-o', '--output', help='Output base dir for the feature parquet files', required=True)

    args = parser.parse_args()

    main(args)
