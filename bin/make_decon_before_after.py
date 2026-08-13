#!/usr/bin/env python

"""QC Tab 4 helper: pair up a "before decon" and "after decon" registration-merged
OME-TIFF (both written by stage_cellprofiler.py, once against the raw image dir and
once against the decon dir - see processes/qc.nf: qc_decon_samples) into one
side-by-side max-projected PNG per channel.

Run once per sampled well by qc_decon_samples. Output filenames
(<plate>_<well>_ch<channel>_before_after.png) are the exact naming convention
tglow.qc.tab_decon.build_before_after_images expects.
"""

import argparse
import glob
import logging
import os

import numpy as np
import tifffile
from matplotlib import pyplot as plt

logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def load_merged_stack(stage_dir):
    """stage_cellprofiler.py --output_format OME_TIFF writes one <field>.ome.tiff per field,
    each a (C, Z, Y, X) stack with every registration-merged channel. QC only needs one
    field, so just take the first one found."""
    paths = sorted(glob.glob(os.path.join(stage_dir, "**", "*.ome.tiff"), recursive=True))
    if not paths:
        raise RuntimeError(f"No .ome.tiff files found under {stage_dir}")

    return tifffile.imread(paths[0])


def max_project_channels(stack):
    """(C, Z, Y, X) -> list of (Y, X) max projections, one per channel."""
    if stack.ndim == 4:
        return [np.max(stack[c], axis=0) for c in range(stack.shape[0])]
    if stack.ndim == 3:
        # Already max-projected upstream (rn_max_project) - (C, Y, X)
        return [stack[c] for c in range(stack.shape[0])]
    raise RuntimeError(f"Unexpected stack shape {stack.shape}, expected (C,Z,Y,X) or (C,Y,X)")


def main(args):
    before_stack = load_merged_stack(args.before)
    after_stack = load_merged_stack(args.after)

    before_channels = max_project_channels(before_stack)
    after_channels = max_project_channels(after_stack)

    n_channels = min(len(before_channels), len(after_channels))
    if len(before_channels) != len(after_channels):
        log.warning(f"Before ({len(before_channels)}) and after ({len(after_channels)}) channel counts differ - showing the common {n_channels}")

    os.makedirs(args.output, exist_ok=True)

    for channel in range(n_channels):
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        for ax, img, title in zip(axes, [before_channels[channel], after_channels[channel]], ["Before decon", "After decon"]):
            vmax = np.percentile(img, 99.5) or 1
            ax.imshow(img, cmap="gray", vmin=0, vmax=vmax)
            ax.set_title(title)
            ax.axis("off")

        fig.suptitle(f"{args.plate} {args.well} - channel {channel}")
        outfile = os.path.join(args.output, f"{args.plate}_{args.well}_ch{channel}_before_after.png")
        fig.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info(f"Wrote {outfile}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build side-by-side before/after decon QC images for one sampled well.")
    parser.add_argument("--before", required=True, help="Dir of registration-merged OME-TIFFs staged from the raw (pre-decon) images")
    parser.add_argument("--after", required=True, help="Dir of registration-merged OME-TIFFs staged from the decon output")
    parser.add_argument("--plate", required=True, help="Reference plate name (used in output filenames)")
    parser.add_argument("--well", required=True, help="Well ID, e.g. A01 (used in output filenames)")
    parser.add_argument("--output", required=True, help="Output directory for the before/after PNGs")
    args = parser.parse_args()
    main(args)
