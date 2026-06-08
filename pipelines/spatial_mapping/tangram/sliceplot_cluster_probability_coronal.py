#!/usr/bin/env python3
"""Create per-slice MERFISH probability scatter plots (coronal view) for Tangram clusters."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sp-h5ad", required=True, help="Spatial AnnData (.h5ad) file")
    parser.add_argument("--metadata-csv", required=True, help="Cell metadata CSV with `brain_section_label`")
    parser.add_argument("--ccf-csv", required=True, help="CCF coordinate CSV with columns x,y,z")
    parser.add_argument("--cluster-csv", required=True,
                        help="Tangram prediction CSV containing per-cluster probabilities per cell")
    parser.add_argument("--output-dir", required=True, help="Directory where slice PDFs will be saved")
    parser.add_argument("--exclude-columns", default="brain_section_label,uniform_density,rna_count_based_density,cluster_id",
                        help="Comma-separated columns to exclude from probability set [default: %(default)s]")
    parser.add_argument("--point-size", type=float, default=2.0,
                        help="Scatter point size [default: %(default)s]")
    parser.add_argument("--padding", type=float, default=0.05,
                        help="Axis padding in mm [default: %(default)s]")
    parser.add_argument("--x-column", default="y_mm",
                        help="Metadata column to use for x-axis (default: %(default)s)")
    parser.add_argument("--y-column", default="z_mm",
                        help="Metadata column to use for y-axis (default: %(default)s)")
    parser.add_argument("--x-label", default="y (Dorsal-Ventral, mm)",
                        help="Label for x-axis (default: %(default)s)")
    parser.add_argument("--y-label", default="z (Medial-Lateral, mm)",
                        help="Label for y-axis (default: %(default)s)")
    return parser.parse_args()


def prepare_anndata(sp_h5ad: Path, metadata_csv: Path, ccf_csv: Path) -> sc.AnnData:
    adata = sc.read_h5ad(sp_h5ad)
    meta = pd.read_csv(metadata_csv, index_col="cell_label")
    ccf = pd.read_csv(ccf_csv, index_col="cell_label")

    adata.obs_names = adata.obs_names.astype(str)
    meta.index = meta.index.astype(str)
    ccf.index = ccf.index.astype(str)

    shared = adata.obs_names.intersection(meta.index).intersection(ccf.index)
    if not len(shared):
        raise ValueError("No overlapping cell barcodes among AnnData, metadata, and CCF tables")

    adata = adata[shared].copy()
    meta = meta.loc[shared]
    ccf = ccf.loc[shared]

    adata.obs["x_mm"] = ccf["x"]
    adata.obs["y_mm"] = ccf["y"]
    adata.obs["z_mm"] = ccf["z"]
    adata.obs["brain_section_label"] = meta["brain_section_label"]
    return adata


def determine_probability_columns(cluster_df: pd.DataFrame, exclude: str) -> list[str]:
    excluded_cols = {col.strip() for col in exclude.split(",") if col.strip()}
    prob_cols = [col for col in cluster_df.columns if col not in excluded_cols]
    if not prob_cols:
        raise ValueError("No probability columns available after applying exclusions")
    return prob_cols


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata_sp = prepare_anndata(Path(args.sp_h5ad), Path(args.metadata_csv), Path(args.ccf_csv))
    cluster_df = pd.read_csv(args.cluster_csv, index_col=0)
    cluster_df.index = cluster_df.index.astype(str)
    cluster_df = cluster_df.reindex(adata_sp.obs_names).fillna(0)

    if args.x_column not in adata_sp.obs.columns or args.y_column not in adata_sp.obs.columns:
        raise ValueError("Specified x/y columns not present in metadata")

    prob_cols = determine_probability_columns(cluster_df, args.exclude_columns)
    global_min = cluster_df[prob_cols].min().min()
    global_max = cluster_df[prob_cols].max().max()

    x_all = adata_sp.obs[args.x_column]
    y_all = adata_sp.obs[args.y_column]
    xmin, xmax = x_all.min() - args.padding, x_all.max() + args.padding
    ymin, ymax = y_all.min() - args.padding, y_all.max() + args.padding

    for cluster in prob_cols:
        cluster_dir = output_dir / f"cluster_{cluster}"
        cluster_dir.mkdir(parents=True, exist_ok=True)
        adata_sp.obs["prob"] = cluster_df[cluster].reindex(adata_sp.obs_names).fillna(0)

        for slice_name in adata_sp.obs["brain_section_label"].unique():
            df = adata_sp.obs[adata_sp.obs["brain_section_label"] == slice_name].copy()
            if df.empty:
                continue
            df = df.sort_values(by="prob", ascending=True)

            colors = np.full(len(df), "#1B1B1B", dtype=object)
            positive = df["prob"] > 0
            if positive.any():
                denom = (global_max - global_min) or 1
                scaled = np.clip((df.loc[positive, "prob"] - global_min) / denom, 0, 1)
                viridis = plt.cm.viridis(scaled)[:, :3]
                colors[positive.to_numpy()] = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in viridis]

            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(df[args.x_column], df[args.y_column], c=colors, s=args.point_size, alpha=1)
            ax.set_title(f"Cluster {cluster} – {slice_name} (n={len(df)})")
            ax.set_xlabel(args.x_label)
            ax.set_ylabel(args.y_label)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect("equal")
            fig.tight_layout()

            sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=global_min, vmax=global_max))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, orientation="vertical")
            cbar.set_label("Cluster Probability", rotation=270, labelpad=15)

            fig.savefig(cluster_dir / f"{slice_name}_cluster_{cluster}.pdf")
            plt.close(fig)


if __name__ == "__main__":
    main()
