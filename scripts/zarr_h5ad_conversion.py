#!/usr/bin/env python3
"""
zarr_h5ad_conversion.py

Convert proseg spatialdata zarr output to files optimized for Seurat v5 import:

  1. proseg-anndata.h5ad  — AnnData with counts, metadata, spatial coords
  2. proseg-seg-vertices.csv.gz — Polygon vertices as (x, y, cell) table,
     ready for CreateSegmentation.data.frame() in R

The vertex table explodes MULTIPOLYGON geometries into ordered vertex rows
grouped by cell name. This bypasses all sf/MULTIPOLYGON issues in R.

Usage:
  python proseg_to_seurat.py <zarr_path> [output_dir]

Example:
  python proseg_to_seurat.py ./proseg-output.zarr ./seurat_ready/

Requirements:
  pip install spatialdata anndata geopandas pandas
"""

import sys
import os
import gzip
import numpy as np
import pandas as pd
import spatialdata


def convert(zarr_path: str, output_dir: str | None = None):
    if output_dir is None:
        output_dir = os.path.dirname(zarr_path) or "."
    os.makedirs(output_dir, exist_ok=True)

    print(f"Reading spatialdata zarr: {zarr_path}")
    sdata = spatialdata.read_zarr(zarr_path)

    # ================================================================
    # 1. Export AnnData as h5ad
    # ================================================================
    adata = sdata.tables["table"]
    h5ad_path = os.path.join(output_dir, "proseg-anndata.h5ad")

    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  obs columns: {list(adata.obs.columns)}")
    print(f"  obsm keys: {list(adata.obsm.keys())}")

    if "spatial" in adata.obsm:
        print(f"  Spatial coords shape: {adata.obsm['spatial'].shape}")

    print(f"Writing h5ad: {h5ad_path}")
    adata.write_h5ad(h5ad_path)

    # ================================================================
    # 2. Export polygon vertices as flat CSV for Seurat
    # ================================================================
    # sdata.shapes["cell_boundaries"] is a GeoDataFrame with cell_id + geometry
    if "cell_boundaries" not in sdata.shapes:
        print("WARNING: No 'cell_boundaries' in sdata.shapes. Skipping vertex export.")
        return h5ad_path

    gdf = sdata.shapes["cell_boundaries"]
    print(f"  Polygons: {len(gdf)} cells")

    # Get cell names from AnnData obs index (matches count matrix)
    cell_names = list(adata.obs_names)

    # Build vertex table: x, y, cell
    # Each polygon's vertices are listed in order, with the cell name repeated
    rows = []
    for idx, (_, row) in enumerate(gdf.iterrows()):
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        # Cell name from AnnData (same row order as shapes)
        if idx < len(cell_names):
            cname = cell_names[idx]
        else:
            cname = str(idx)

        # Handle both Polygon and MultiPolygon
        if geom.geom_type == "MultiPolygon":
            # Take exterior of all sub-polygons
            for poly in geom.geoms:
                coords = np.array(poly.exterior.coords)
                for x, y in coords:
                    rows.append((x, y, cname))
        elif geom.geom_type == "Polygon":
            coords = np.array(geom.exterior.coords)
            for x, y in coords:
                rows.append((x, y, cname))

    vtx_df = pd.DataFrame(rows, columns=["x", "y", "cell"])
    vtx_path = os.path.join(output_dir, "proseg-seg-vertices.csv.gz")

    print(f"  Vertices: {len(vtx_df)} rows")
    print(f"Writing vertex table: {vtx_path}")
    vtx_df.to_csv(vtx_path, index=False, compression="gzip")

    # ================================================================
    # 3. Export centroids as CSV (backup, also in h5ad obsm)
    # ================================================================
    if "spatial" in adata.obsm:
        cent_df = pd.DataFrame(
            adata.obsm["spatial"],
            columns=["x", "y"],
            index=adata.obs_names,
        )
        cent_df["cell"] = cent_df.index
        cent_path = os.path.join(output_dir, "proseg-centroids.csv.gz")
        print(f"Writing centroids: {cent_path}")
        cent_df.to_csv(cent_path, index=False, compression="gzip")

    print("Done.")
    return h5ad_path


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python proseg_to_seurat.py <zarr_path> [output_dir]")
        sys.exit(1)
    
    zarr_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    convert(zarr_path, output_dir)
