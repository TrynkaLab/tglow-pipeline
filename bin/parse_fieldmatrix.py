#!/usr/bin/env python

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tglow.io.image_query import ImageQuery


def make_field_matrix(df, drop_z=True):
    """
    df: DataFrame with columns ['field', 'x', 'y', 'z']
    drop_z: if True -> return 2D (y, x) matrix
            if False -> return 3D (z, y, x) matrix
    Returns (mat, df_with_indices)
    """
    # Work on a copy
    df = df.copy()
    df = df.sort_values(['y', 'x', 'z']).drop_duplicates(subset=['x', 'y', 'z', 'field'])

    xs = np.sort(df['x'].unique())
    ys = np.flip(np.sort(df['y'].unique()))

    x_to_idx = {v: i for i, v in enumerate(xs)}
    y_to_idx = {v: i for i, v in enumerate(ys)}

    # Add index columns to df (i = row (y), j = col (x))
    df['index_y'] = df['y'].map(y_to_idx)
    df['index_x'] = df['x'].map(x_to_idx)

    if drop_z:
        # 2D matrix: rows = y (i), cols = x (j)
        mat = np.full((len(ys), len(xs)), np.nan)
        for _, row in df.iterrows():
            i = int(row['index_y'])
            j = int(row['index_x'])
            mat[i, j] = row['field']
    else:
        # 3D matrix: axes = (z, y, x)
        zs = np.sort(df['z'].unique())
        z_to_idx = {v: k for k, v in enumerate(zs)}
        df['index_z'] = df['z'].map(z_to_idx)
        mat = np.full((len(zs), len(ys), len(xs)), np.nan)
        for _, row in df.iterrows():
            k = int(row['index_z'])
            i = int(row['index_y'])
            j = int(row['index_x'])
            mat[k, i, j] = row['field']

    return mat, df

def main(args):
    input_file = Path(args.file)
    drop_z = args.drop_z
    outdir = Path(args.output_prefix)

    if not input_file.exists():
        print(f"ERROR: input file not found: {input_file}", file=sys.stderr)
        sys.exit(2)

    outdir.mkdir(parents=True, exist_ok=True)

    with input_file.open() as fh:
        js = json.load(fh)

    all_indices = []
    saved_matrix = False
    for i_well, well in enumerate(js.get('wells', [])):
        images = well.get('images', [])
        if len(images) == 0:
            continue

        df = pd.DataFrame(images)
        pos = pd.json_normalize(df['position'])
        coords = pd.DataFrame({
            'field': df['field'],
            'x': pos.get('x.value'),
            'y': pos.get('y.value'),
            'z': pos.get('z.value'),
        })

        if drop_z:
            coords['z'] = 0

        coords_unique = coords.drop_duplicates()

        iq = ImageQuery(js['plate']['name'], well['row'], well['col'], None)

        mat, indices = make_field_matrix(coords_unique, drop_z=False)
        indices['well'] = iq.get_well_id()
        indices['plate'] = iq.plate
        indices['row'] = iq.get_row_letter()
        indices['col'] = iq.col

        #indices['index_y'] = indices['index_y']
        #indices['index_x'] = indices['index_x']
        #indices['index_z'] = indices['index_z']

        # Save example matrix for first processed well
        if not saved_matrix:
            if mat.ndim == 3:
                np.savetxt(outdir / 'xy_field_matrix_well1_z1.tsv', mat[0, :, :], delimiter='\t', fmt='%g')
            else:
                np.savetxt(outdir / 'xy_field_matrix_well1_z1.tsv', mat, delimiter='\t', fmt='%g')
            saved_matrix = True

        all_indices.append(indices)
        

    if all_indices:
        all_indices_df = pd.concat(all_indices, ignore_index=True)
        all_indices_df.to_csv(outdir / 'field_indices_per_well.tsv', sep='\t', index=False)
        print(f"Wrote {len(all_indices_df)} index rows to {outdir / 'field_indices_per_well.tsv'}")
    else:
        print("No indices collected, no output written.")


if __name__ == "__main__":
    
    p = argparse.ArgumentParser(description="Parse Index.json into field matrices and indices TSV")
    p.add_argument("--file", "-f", required=True, help="Path to Index.json")
    p.add_argument("--drop_z", action="store_true", help="Collapse z dimension (produce 2D matrices)")
    p.add_argument("--output_prefix", "-o", default=".", help="Output directory (will be created if needed)")
    args = p.parse_args()
    
    main(args)