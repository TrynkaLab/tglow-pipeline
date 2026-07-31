#!/usr/bin/env python

"""Calculate per-plate/per-channel scaling factors and sigmoid soft-threshold
parameters from cell- and image-level feature h5ads.

Python translation of scaling_script.r. See scaling_script.md for the input
format this expects (h5ads instead of a tglow object, channel map and control
list as input files instead of being baked into the tglow object).
"""

import argparse
import logging

import numpy as np
import pandas as pd
import anndata as ad

from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_utils import sigmoid_params

logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

REQUIRED_CHANNEL_MAP_COLS = [
    "channel",
    "name",
    "rep_features_dynamic",
    "rep_features_offset",
    "rep_features_scale_lower",
    "rep_features_scale_upper",
    "control_population",
]

# Optional channel_map columns; missing entirely or blank per-row both default to False.
TRUE_VALUES = {"true", "t", "1", "yes", "y"}


def _to_bool(val):
    """Best-effort coercion of a channel_map cell to bool, defaulting missing/NaN to False."""
    if pd.isna(val):
        return False
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    if isinstance(val, (int, float, np.integer, np.floating)):
        return bool(val)
    return str(val).strip().lower() in TRUE_VALUES


def load_channel_map(path):
    """Load the channel map TSV and drop any ref_plate/plate columns.

    The channel map describes one reference plate/cycle's channel layout
    (effectively channel_indices.tsv), extended with the rep_features_*/
    control_population columns. It gets broadcast onto every real plate found
    in the data by build_scaling_index, so any ref_plate/plate columns it
    carries in from channel_indices.tsv are dropped here (R:
    channel.map.orig[,-c(1,2)]).

    Three optional per-channel boolean columns are also recognized (missing
    entirely, or blank per row, both default to False):
    - skip_scaling: force scale_factor=1 for this channel regardless of what
      rep_features_dynamic/rep_features_offset would otherwise compute.
    - skip_sigmoid: always use the "no soft-threshold" sigmoid fallback for
      this channel, regardless of whether rep_features_scale_lower/upper are
      set (e.g. brightfield, where dynamic-range scaling is still measured
      but a sigmoid threshold isn't meaningful).
    - skip_plate_offsets: ignore the control-derived plate_offset for this
      channel (forced to 1 for every plate) -- the channel is still scaled to
      its measured dynamic range, just without cross-plate relative
      adjustment. plate_feature_mean/plate_ncells are still computed/reported.
    """
    channel_map = pd.read_csv(path, sep="\t")

    missing = [c for c in REQUIRED_CHANNEL_MAP_COLS if c not in channel_map.columns]
    if missing:
        raise ValueError(f"channel_map is missing required columns: {missing}")

    channel_map = channel_map.drop(columns=["ref_plate", "plate"], errors="ignore")
    return channel_map


def get_feature_series(adata, feature):
    """Extract a single feature column from adata.X as a dense 1D array."""
    idx = adata.var_names.get_loc(feature)
    col = adata.X[:, idx]

    if hasattr(col, "toarray"):
        col = col.toarray()

    return np.asarray(col).ravel()


def build_well_column(obs, row_col, col_col):
    """Build 'A01'-style well ids from row/col columns, matching ImageQuery.get_well_id()."""
    row_letter = obs[row_col].astype(int).astype(str).map(ImageQuery.ID_TO_ROW)
    col_str = obs[col_col].astype(int).astype(str).str.zfill(2)
    return row_letter + col_str


def build_scaling_index(channel_map, plates):
    """Cross join channel_map x plates, indexed by '<plate>-<channel>'."""
    scaling_index = pd.concat(
        [channel_map.assign(ref_plate=plate) for plate in plates],
        ignore_index=True,
    )
    scaling_index.index = scaling_index["ref_plate"].astype(str) + "-" + scaling_index["channel"].astype(str)
    return scaling_index


def compute_dynamic_range(scaling_index, channel_map, cell_adata, plate_col, q2):
    """Per channel: total and per-plate quantile(q2) of the dynamic-range feature (R lines 147-168)."""
    plates_series = cell_adata.obs[plate_col].astype(str)

    for _, ch_row in channel_map.iterrows():
        channel = ch_row["channel"]
        feature = ch_row["rep_features_dynamic"]
        mask = scaling_index["channel"] == channel

        if pd.isna(feature):
            scaling_index.loc[mask, "max_scale_total"] = 1.0
            scaling_index.loc[mask, "max_scale_plate"] = 1.0
            continue

        values = get_feature_series(cell_adata, feature)
        df = pd.DataFrame({"plate": plates_series.values, feature: values})

        scaling_index.loc[mask, "max_scale_total"] = df[feature].quantile(q2)

        per_plate_q = df.groupby("plate")[feature].quantile(q2)
        for plate in scaling_index.loc[mask, "ref_plate"].unique():
            if plate in per_plate_q.index:
                key = f"{plate}-{channel}"
                scaling_index.loc[key, "max_scale_plate"] = per_plate_q[plate]

    return scaling_index


def load_controls(path):
    """Load the <plate> <well> <control_type> controls TSV."""
    controls = pd.read_csv(path, sep="\t", dtype={"plate": str, "well": str})
    required = ["plate", "well", "control_type"]
    missing = [c for c in required if c not in controls.columns]
    if missing:
        raise ValueError(f"controls file is missing required columns: {missing}")
    return controls


def _merge_controls(obs, controls, plate_col, row_col, col_col):
    """Inner-merge controls onto obs via plate+well, returning only matching (control) rows.

    Preserves obs's original row index (merge() would otherwise reset it),
    since callers use this index to select back into the source AnnData's
    feature arrays.

    control_population matching against control_type is intentionally not
    implemented yet (see scaling_script.md) -- every row present in the
    controls file counts as a control, regardless of its control_type value.
    """
    df = obs.copy()
    df["plate"] = df[plate_col].astype(str)
    df["well"] = build_well_column(df, row_col, col_col)
    df["_obs_index"] = df.index

    merged = df.merge(controls, on=["plate", "well"], how="inner", suffixes=(None, "_control"))
    merged = merged.set_index("_obs_index")
    return merged


def apply_skip_plate_offsets(scaling_index, channel_map):
    """Force plate_offset=1 for channels flagged skip_plate_offsets=True.

    The channel still gets scaled to its dynamic range (rep_features_dynamic),
    just without the cross-plate relative adjustment derived from control means.
    """
    if "skip_plate_offsets" not in channel_map.columns:
        return scaling_index

    for _, ch_row in channel_map.iterrows():
        if _to_bool(ch_row.get("skip_plate_offsets")):
            mask = scaling_index["channel"] == ch_row["channel"]
            scaling_index.loc[mask, "plate_offset"] = 1.0

    return scaling_index


def compute_plate_offsets(scaling_index, channel_map, cell_adata, controls_merged, plate_col):
    """Per channel: plate offsets from control-cell means (R lines 173-224)."""
    for _, ch_row in channel_map.iterrows():
        channel = ch_row["channel"]
        feature = ch_row["rep_features_offset"]
        mask = scaling_index["channel"] == channel

        if pd.isna(feature):
            scaling_index.loc[mask, "plate_feature_mean"] = np.nan
            scaling_index.loc[mask, "plate_offset"] = 1.0
            continue

        values = get_feature_series(cell_adata, feature)
        cur_df = pd.DataFrame({
            "plate": cell_adata.obs[plate_col].astype(str).values,
            feature: values,
        }, index=cell_adata.obs.index).loc[controls_merged.index]

        ncells = cur_df.dropna(subset=[feature]).groupby("plate").size()
        plate_means = cur_df.groupby("plate")[feature].mean()

        if plate_means.empty:
            continue

        plate_offset = plate_means / plate_means.min()

        for plate, mean_val in plate_means.items():
            key = f"{plate}-{channel}"
            scaling_index.loc[key, "plate_feature_mean"] = mean_val
            scaling_index.loc[key, "plate_offset"] = plate_offset[plate]
            scaling_index.loc[key, "plate_ncells"] = ncells.get(plate, 0)

    # If plate_offset is NA (no control cells for that plate/channel), default to 1
    scaling_index["plate_offset"] = scaling_index["plate_offset"].fillna(1.0)

    # skip_plate_offsets: keep plate_feature_mean/plate_ncells (still measured/reported),
    # but ignore the offset when deriving the scale factor below -- each plate's scale
    # factor then reflects only its own dynamic range, with no cross-plate adjustment.
    scaling_index = apply_skip_plate_offsets(scaling_index, channel_map)

    # Normalize the base plate scale
    scaling_index["max_scale_plate_norm"] = scaling_index["max_scale_plate"] / scaling_index["plate_offset"]

    # Base scale: highest normalized saturation point per channel
    for channel in scaling_index["channel"].unique():
        mask = scaling_index["channel"] == channel
        scaling_index.loc[mask, "base_scale"] = scaling_index.loc[mask, "max_scale_plate_norm"].max()

    scaling_index["scale_factor"] = scaling_index["base_scale"] * scaling_index["plate_offset"]
    scaling_index["scale_factor"] = scaling_index["scale_factor"].fillna(1.0)

    return scaling_index


def apply_skip_scaling(scaling_index, channel_map):
    """Force scale_factor=1 (no scaling applied) for channels flagged skip_scaling=True.

    Overrides whatever was computed from rep_features_dynamic/rep_features_offset,
    for channels (e.g. brightfield) that get measured like any other channel but
    should not have a scale factor applied.
    """
    if "skip_scaling" not in channel_map.columns:
        return scaling_index

    for _, ch_row in channel_map.iterrows():
        if _to_bool(ch_row.get("skip_scaling")):
            mask = scaling_index["channel"] == ch_row["channel"]
            scaling_index.loc[mask, "scale_factor"] = 1.0

    return scaling_index


def compute_sigmoid_params(scaling_index, channel_map, image_adata, controls_merged, plate_col, sigmoid_tol, uint_max):
    """Per channel: fit sigmoid bias/slope from control-image lower/upper quantiles (R lines 236-292)."""
    for _, ch_row in channel_map.iterrows():
        channel = ch_row["channel"]
        lower_feature = ch_row["rep_features_scale_lower"]
        upper_feature = ch_row["rep_features_scale_upper"]
        mask = scaling_index["channel"] == channel

        if pd.isna(lower_feature) or _to_bool(ch_row.get("skip_sigmoid")):
            # No feature set, or sigmoid explicitly disabled (e.g. brightfield):
            # fall back to a curve that always evaluates ~1
            par = sigmoid_params(-99, 1, tol=1e-16)
            scaling_index.loc[mask, "sigmoid_bias"] = par["bias"]
            scaling_index.loc[mask, "sigmoid_slope"] = par["slope"]
            scaling_index.loc[mask, "sigmoid_tol"] = par["tol"]
            scaling_index.loc[mask, "sigmoid_x1"] = -99
            scaling_index.loc[mask, "sigmoid_x2"] = 1
            continue

        lower_values = get_feature_series(image_adata, lower_feature)
        upper_values = get_feature_series(image_adata, upper_feature)

        cur_df = pd.DataFrame({
            "plate": image_adata.obs[plate_col].astype(str).values,
            "lower": lower_values,
            "upper": upper_values,
        }, index=image_adata.obs.index).loc[controls_merged.index]

        lower_q = cur_df.groupby("plate")["lower"].quantile(0.95) * uint_max
        upper_q = cur_df.groupby("plate")["upper"].median() * uint_max

        for plate in lower_q.index:
            x1 = lower_q[plate]
            x2 = upper_q.get(plate, np.nan)
            if pd.isna(x1) or pd.isna(x2):
                continue

            par = sigmoid_params(x1, x2, tol=sigmoid_tol)
            key = f"{plate}-{channel}"
            scaling_index.loc[key, "sigmoid_bias"] = par["bias"]
            scaling_index.loc[key, "sigmoid_slope"] = par["slope"]
            scaling_index.loc[key, "sigmoid_tol"] = par["tol"]
            scaling_index.loc[key, "sigmoid_x1"] = x1
            scaling_index.loc[key, "sigmoid_x2"] = x2

    return scaling_index


def fill_missing_sigmoid(scaling_index):
    """Hotfix: fill missing sigmoid slope/bias with the channel's mean (R lines 294-307)."""
    missing_mask = scaling_index["sigmoid_slope"].isna()

    if missing_mask.any():
        log.warning(
            "Not all plates have a scaling index, this can happen if plates miss control "
            "samples. As a hack, setting to the mean of the plates"
        )
        for channel in scaling_index.loc[missing_mask, "channel"].unique():
            channel_mask = scaling_index["channel"] == channel
            selector = channel_mask & scaling_index["sigmoid_slope"].isna()

            scaling_index.loc[selector, "sigmoid_slope"] = scaling_index.loc[channel_mask, "sigmoid_slope"].mean()
            scaling_index.loc[selector, "sigmoid_bias"] = scaling_index.loc[channel_mask, "sigmoid_bias"].mean()

    return scaling_index


def compute_channel_medians(scaling_index):
    """Per-channel median sigmoid slope/bias, broadcast to all plates of that channel (R lines 309-318)."""
    for channel in scaling_index["channel"].unique():
        mask = scaling_index["channel"] == channel
        scaling_index.loc[mask, "sigmoid_slope_mean"] = scaling_index.loc[mask, "sigmoid_slope"].median()
        scaling_index.loc[mask, "sigmoid_bias_mean"] = scaling_index.loc[mask, "sigmoid_bias"].median()

    return scaling_index


def _write_channel_value_file(scaling_index, column, path):
    """Write '<plate>_ch<channel>=<value>' tokens, space separated on one line (R lines 322-334)."""
    tokens = (
        scaling_index["ref_plate"].astype(str)
        + "_ch"
        + scaling_index["channel"].astype(str)
        + "="
        + scaling_index[column].astype(str)
    )
    with open(path, "w") as fh:
        fh.write(" ".join(tokens))


def write_outputs(scaling_index, output_dir):
    scaling_index.to_csv(f"{output_dir}/scaling_index.tsv", sep="\t", index=True, index_label="scaling_key")
    _write_channel_value_file(scaling_index, "scale_factor", f"{output_dir}/scaling_factors.txt")
    _write_channel_value_file(scaling_index, "sigmoid_bias", f"{output_dir}/sigmoid_bias.txt")
    _write_channel_value_file(scaling_index, "sigmoid_slope", f"{output_dir}/sigmoid_slope.txt")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate per-plate/per-channel scaling factors and sigmoid parameters "
                    "from cell- and image-level feature h5ads."
    )
    parser.add_argument("--cell_h5ad", required=True, help="Path to cell-level features h5ad")
    parser.add_argument("--image_h5ad", required=True, help="Path to image-level features h5ad")
    parser.add_argument("--channel_map", required=True, help="TSV describing the channel layout, extended with rep_features_*/control_population columns")
    parser.add_argument("--controls", required=True, help="TSV with <plate> <well> <control_type> columns")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--q2", type=float, default=0.99, help="Quantile (0-1) used for dynamic range")
    parser.add_argument("--sigmoid_tol", type=float, default=1e-3, help="Sigmoid tolerance for scale_lower/upper fit")
    parser.add_argument("--uint_max", type=int, default=65535, help="Max value of the raw intensity range")
    parser.add_argument("--plate_col", default="plate")
    parser.add_argument("--row_col", default="row")
    parser.add_argument("--col_col", default="col")
    parser.add_argument("--field_col", default="field")
    parser.add_argument("--cell_col", default="cell_number")
    args = parser.parse_args()

    channel_map = load_channel_map(args.channel_map)
    controls = load_controls(args.controls)

    cell_adata = ad.read_h5ad(args.cell_h5ad)
    image_adata = ad.read_h5ad(args.image_h5ad)

    plates = sorted(set(cell_adata.obs[args.plate_col].astype(str)) | set(image_adata.obs[args.plate_col].astype(str)))
    log.info(f"Found {len(plates)} plates: {plates}")

    scaling_index = build_scaling_index(channel_map, plates)

    scaling_index = compute_dynamic_range(scaling_index, channel_map, cell_adata, args.plate_col, args.q2)

    cell_controls = _merge_controls(cell_adata.obs, controls, args.plate_col, args.row_col, args.col_col)
    scaling_index = compute_plate_offsets(scaling_index, channel_map, cell_adata, cell_controls, args.plate_col)
    scaling_index = apply_skip_scaling(scaling_index, channel_map)

    image_controls = _merge_controls(image_adata.obs, controls, args.plate_col, args.row_col, args.col_col)
    scaling_index = compute_sigmoid_params(scaling_index, channel_map, image_adata, image_controls, args.plate_col, args.sigmoid_tol, args.uint_max)
    scaling_index = fill_missing_sigmoid(scaling_index)
    scaling_index = compute_channel_medians(scaling_index)

    write_outputs(scaling_index, args.output)
    log.info("Done")


if __name__ == "__main__":
    main()
