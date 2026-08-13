#!/usr/bin/env python

"""Build the self-contained QC report HTML page (Tab 1 always, tabs 2-6 optional -
see subworkflows/qc_report.nf for which artifacts get passed for which tabs). Thin
CLI wrapper - all the actual work is in tglow.qc.render.build_report.

Nextflow passes disabled optional inputs as placeholder paths (e.g. file("NO_X")),
following the same convention used throughout the rest of the pipeline - a
placeholder simply won't exist on disk, so --resolve_optional treats it as absent.
"""

import argparse
import json
import logging
import os

from tglow.qc.render import build_report

logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def resolve_optional(path):
    """None or a Nextflow NO_X placeholder (doesn't exist on disk) -> None."""
    if path is None or not os.path.exists(path):
        return None
    return path


def main(args):
    ff_params = json.loads(args.ff_params) if args.ff_params else {}

    build_report(
        output_path=args.output,
        measurements_dir=args.measurements_dir,
        manifest_path=args.manifest,
        blacklist_path=resolve_optional(args.blacklist),
        registration_manifest_path=resolve_optional(args.registration_manifest),
        registration_images_dir=resolve_optional(args.registration_images_dir),
        qc_regcor=args.qc_regcor,
        qc_registration_pattern=args.qc_registration_pattern,
        qc_plate_format=args.qc_plate_format,
        qc_n_sample_registration=args.qc_n_sample_registration,
        show_flatfield=args.show_flatfield,
        flatfields_dir=resolve_optional(args.flatfields_dir),
        ff_global_flatfield=args.ff_global_flatfield,
        ff_params=ff_params,
        show_decon=args.show_decon,
        decon_samples_dir=resolve_optional(args.decon_samples_dir),
        show_scaling=args.show_scaling,
        scaling_index_path=resolve_optional(args.scaling_index),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the tglow QC report HTML page.")

    parser.add_argument("--measurements_dir", required=True, help="Per-plate dir of measure_intensity parquet output (as staged by stage_as_plate)")
    parser.add_argument("--manifest", required=True, help="rn_manifest.tsv path")
    parser.add_argument("--output", required=True, help="Output HTML path")

    parser.add_argument("--blacklist", default=None, help="rn_blacklist.tsv path, or a NO_X placeholder if not set")
    parser.add_argument("--qc_regcor", type=float, default=0.6)
    parser.add_argument("--qc_registration_pattern", default="registration_corr")
    parser.add_argument("--qc_plate_format", default="auto")
    parser.add_argument("--qc_n_sample_registration", type=int, default=10)

    # Tab 2: registration
    parser.add_argument("--registration_manifest", default=None, help="rn_manifest_registration.tsv path, or a NO_X placeholder")
    parser.add_argument("--registration_images_dir", default=None, help="Root dir of registration PNGs (rr../registration), or a NO_X placeholder")

    # Tab 3: flatfields
    parser.add_argument("--show_flatfield", action="store_true", default=False)
    parser.add_argument("--flatfields_dir", default=None, help="Root dir of flatfield PNGs")
    parser.add_argument("--ff_global_flatfield", action="store_true", default=False)
    parser.add_argument("--ff_params", default=None, help="JSON dict of ff_* param values, for the Tab 3 summary table")

    # Tab 4: decon
    parser.add_argument("--show_decon", action="store_true", default=False)
    parser.add_argument("--decon_samples_dir", default=None, help="Output dir of the qc_decon_samples process (before/after PNGs)")

    # Tab 6: scaling factors
    parser.add_argument("--show_scaling", action="store_true", default=False)
    parser.add_argument("--scaling_index", default=None, help="scaling_index.tsv path")

    args = parser.parse_args()
    main(args)
