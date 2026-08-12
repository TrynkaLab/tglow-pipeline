#!/usr/bin/env python

"""Assign a short, deterministic index (P1, P2, ...) to every plate in the core manifest.

Indices are assigned by sorted plate name, not manifest row order or arrival
order, so the mapping is stable regardless of how Nextflow schedules anything
downstream - it depends only on the fixed set of plate names in the manifest.
"""

import argparse
import logging

import pandas as pd

logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def build_plate_index(manifest_path):
    """Return a plate -> plate_id (P1, P2, ...) dict, sorted by plate name."""
    manifest = pd.read_csv(manifest_path, sep="\t")
    plates = sorted(manifest["plate"].unique())
    return {plate: f"P{i}" for i, plate in enumerate(plates, start=1)}


def main(args):
    plate_index = build_plate_index(args.manifest)

    with open(args.output, "w") as f:
        f.write("plate\tplate_id\n")
        for plate, plate_id in plate_index.items():
            f.write(f"{plate}\t{plate_id}\n")

    log.info(f"Wrote {len(plate_index)} plate indices to {args.output}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Assign a short, deterministic index to every plate in the core manifest.")

    parser.add_argument('--manifest', help='Path to the core manifest TSV (must have a "plate" column)', required=True)
    parser.add_argument('--output', help='Output path for the plate\\tplate_id TSV', required=True)

    args = parser.parse_args()

    main(args)
