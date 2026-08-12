#!/usr/bin/env python

"""Aggregate per-well CellProfiler output into per-plate cell-level and image-level parquet files.

Python translation of read_cellprofiler_fileset_b/add_global_ids in
tglow-r/R/io_cellprofiler.r, adapted to run per-plate instead of once over a
whole experiment - CellProfiler's own ImageNumber/ObjectNumber only need to be
unique within one well's output, so this rebuilds a global id per plate from
a short plate index (--plate_id, e.g. "P1") and a well index derived locally
by sorting the plate's well zips alphabetically (both deterministic
regardless of Nextflow's actual, non-deterministic task-completion order).

Each well's CellProfiler output is expected as a self-contained
<plate>_<well>.zip (as written by processes/cellprofiler.nf's `cellprofiler`/
`finalize_and_cellprofiler`), containing one _cell.txt (parent objects), one
_Image.txt (image-level data), optionally _Experiment.txt/_Object
relationships.txt (ignored), and one file per child object type whose rows
are matched to their parent cell via --parent_col and merged onto it.
"""

import argparse
import glob
import logging
import os
import re
import zipfile

import pandas as pd

logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def find_well_zips(input_dir, plate):
    """Every <plate>_<well>.zip in input_dir, sorted by filename (== sorted by well)."""
    return sorted(glob.glob(os.path.join(input_dir, f"{plate}_*.zip")))


def parse_well_from_zipname(zip_path, plate):
    name = os.path.basename(zip_path)
    prefix = f"{plate}_"
    if not name.startswith(prefix) or not name.endswith(".zip"):
        raise ValueError(f"Zip filename '{name}' does not match the expected '{plate}_<well>.zip' pattern")
    return name[len(prefix):-len(".zip")]


def classify_member(name, child_pattern):
    """Classify a zip member as ('cell'|'image'|'skip', None) or ('child', object_name)."""
    base = os.path.basename(name)
    if base.endswith("_cell.txt"):
        return "cell", None
    if base.endswith("_Image.txt"):
        return "image", None
    if base.endswith("_Experiment.txt"):
        return "skip", None
    if "Object relationships" in base:
        return "skip", None
    match = re.match(child_pattern, base)
    if match:
        return "child", match.group(1)
    return "skip", None


def read_member(zf, name):
    with zf.open(name) as fh:
        try:
            return pd.read_csv(fh, sep="\t")
        except pd.errors.EmptyDataError:
            return pd.DataFrame()


def add_global_ids(df, fileset_id):
    """Add <col>_Global columns for every ImageNumber/ObjectNumber/Object_Number/Parent
    column, mirroring add_global_ids() in tglow-r's io_cellprofiler.r."""
    df = df.copy()

    image_cols = [c for c in df.columns if "ImageNumber" in c]
    for col in image_cols:
        df[f"{col}_Global"] = df[col].map(lambda v: f"FS{fileset_id}_I{v}" if pd.notna(v) else pd.NA)

    if image_cols:
        img_values = df[image_cols[0]]
        object_cols = [c for c in df.columns if ("ObjectNumber" in c) or ("Object_Number" in c) or ("Parent" in c)]
        for col in object_cols:
            df[f"{col}_Global"] = [
                f"FS{fileset_id}_I{img}_O{obj}" if pd.notna(img) and pd.notna(obj) else pd.NA
                for img, obj in zip(img_values, df[col])
            ]

    return df


def merge_child_object(cells, cur, obj, parent_col, merge_strategy, zip_path):
    """Aggregate one child object's rows onto `cells` (which already carries a
    "_merge_key" == cell_ImageNumber_cell_ObjectNumber column), returning the
    updated `cells`."""

    if len(cur.columns) == 0:
        log.warning(f"{zip_path}: child object '{obj}' file has no header, skipping (no columns to add)")
        cells[f"{obj}_QC_Object_Count"] = 0
        return cells

    if parent_col not in cur.columns:
        raise ValueError(
            f"Parent col '{parent_col}' not found for child object '{obj}' in {zip_path}. "
            f"Update --parent_col or the CellProfiler pipeline."
        )

    # Unlike R's `cur[[parent.col]] != 0` (which, via R's NA-preserving logical
    # indexing, would keep a row of NA for unassigned/missing children), also
    # drop NaN parent values outright - there's nothing meaningful to merge.
    cur = cur[(cur[parent_col] != 0) & cur[parent_col].notna()].copy()

    if len(cur) == 0:
        log.warning(f"{zip_path}: child object '{obj}' has no rows assigned to a parent, filling NA")
        for col in cur.columns:
            cells[f"{obj}_{col}"] = pd.NA
        cells[f"{obj}_QC_Object_Count"] = 0
        return cells

    cur.columns = [f"{obj}_{c}" for c in cur.columns]
    child_key = cur[f"{obj}_ImageNumber"].astype(str) + "_" + cur[f"{obj}_{parent_col}"].astype(str)
    child_key.name = "_merge_key"

    # Drop per-child object-identifier columns before aggregating - meaningless
    # once merged onto the parent (matches R's read_cellprofiler_fileset_b exclude.cols).
    exclude_cols = {c for c in cur.columns if "ObjectNumber" in c or "Number_Object_Number" in c}
    agg_cols = [c for c in cur.columns if c not in exclude_cols]
    numeric_cols = [c for c in agg_cols if pd.api.types.is_numeric_dtype(cur[c])]
    categorical_cols = [c for c in agg_cols if c not in numeric_cols]

    grouped = cur.groupby(child_key)
    parts = []
    if numeric_cols:
        parts.append(grouped[numeric_cols].agg(merge_strategy))
    if categorical_cols:
        # Categorical columns aren't averageable - keep the value if every row
        # in the group agrees, otherwise NA (generalizes R, which only ever
        # applied 'mean' and implicitly assumed every column was numeric).
        parts.append(grouped[categorical_cols].agg(lambda s: s.iloc[0] if s.nunique(dropna=False) == 1 else pd.NA))

    agg = pd.concat(parts, axis=1) if parts else pd.DataFrame(index=grouped.size().index)
    agg[f"{obj}_QC_Object_Count"] = grouped.size()
    agg.index.name = "_merge_key"

    return cells.merge(agg.reset_index(), how="left", on="_merge_key")


def process_well(zip_path, plate, well, well_idx, plate_id, parent_col, child_pattern, merge_strategy):
    """Return (cells_df, image_df) for one well, or (None, None) if it has nothing usable."""

    fileset_id = f"{plate_id}W{well_idx}"

    with zipfile.ZipFile(zip_path) as zf:
        cell_member = image_member = None
        child_members = {}

        for member in zf.namelist():
            kind, obj = classify_member(member, child_pattern)
            if kind == "cell":
                cell_member = member
            elif kind == "image":
                image_member = member
            elif kind == "child":
                if obj in child_members:
                    log.warning(f"{zip_path}: multiple files matched child object '{obj}', using {member} (overwriting {child_members[obj]})")
                child_members[obj] = member

        if cell_member is None or image_member is None:
            log.warning(f"{zip_path}: missing _cell.txt or _Image.txt, skipping well")
            return None, None

        cells = read_member(zf, cell_member)
        image = read_member(zf, image_member)

        if len(cells) == 0:
            log.warning(f"{zip_path}: _cell.txt has no rows, skipping well")
            return None, None

        cells.columns = [f"cell_{c}" for c in cells.columns]
        cells["_merge_key"] = cells["cell_ImageNumber"].astype(str) + "_" + cells["cell_ObjectNumber"].astype(str)

        for obj, member in sorted(child_members.items()):
            cur = read_member(zf, member)
            cells = merge_child_object(cells, cur, obj, parent_col, merge_strategy, zip_path)

    cells = cells.drop(columns=["_merge_key"])

    cells = add_global_ids(cells, fileset_id)
    image = add_global_ids(image, fileset_id)

    for df in (cells, image):
        df["plate"] = plate
        df["well"] = well
        df["well_idx"] = well_idx

    cells["global_cell_id"] = cells.get("cell_ObjectNumber_Global")
    image["global_image_id"] = image.get("ImageNumber_Global")

    return cells, image


def main(args):
    zip_paths = find_well_zips(args.input, args.plate)
    if not zip_paths:
        raise RuntimeError(f"No '{args.plate}_*.zip' files found in {args.input}")

    cells_parts = []
    image_parts = []

    for well_idx, zip_path in enumerate(zip_paths, start=1):
        well = parse_well_from_zipname(zip_path, args.plate)
        log.info(f"Processing {args.plate} well {well} ({well_idx}/{len(zip_paths)}): {zip_path}")

        cells, image = process_well(
            zip_path, args.plate, well, well_idx, args.plate_id,
            args.parent_col, args.child_pattern, args.merge_strategy
        )
        if cells is not None:
            cells_parts.append(cells)
        if image is not None:
            image_parts.append(image)

    if not cells_parts:
        log.warning(f"No wells produced any cell rows for plate {args.plate}")

    all_cells = pd.concat(cells_parts, ignore_index=True) if cells_parts else pd.DataFrame()
    all_image = pd.concat(image_parts, ignore_index=True) if image_parts else pd.DataFrame()

    all_cells.to_parquet(args.output_cells, index=False)
    all_image.to_parquet(args.output_image, index=False)

    log.info(f"Wrote {args.output_cells} ({len(all_cells)} rows, {len(all_cells.columns)} cols)")
    log.info(f"Wrote {args.output_image} ({len(all_image)} rows, {len(all_image.columns)} cols)")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Aggregate per-well CellProfiler output into per-plate cell-level and image-level parquet files.")

    parser.add_argument('--input', required=True, help="Directory containing this plate's <plate>_<well>.zip files")
    parser.add_argument('--plate', required=True, help='Plate to process')
    parser.add_argument('--plate_id', required=True, help='Short deterministic plate index (e.g. P1), used to build globally unique cell/image ids')
    parser.add_argument('--merge_strategy', default='mean', choices=['mean', 'median', 'sum'], help='How to aggregate numeric child-object columns onto their parent cell')
    parser.add_argument('--parent_col', default='Parent_cell', help='Column used to match a child object row to its parent cell')
    parser.add_argument('--child_pattern', default=r'^.*_([a-zA-Z]+\d*)\.txt$', help='Regex (first capture group = object name) used to identify a file as a child object and name it')
    parser.add_argument('--output_cells', required=True, help='Output parquet path for the cell-level table')
    parser.add_argument('--output_image', required=True, help='Output parquet path for the image-level table')

    args = parser.parse_args()

    main(args)
