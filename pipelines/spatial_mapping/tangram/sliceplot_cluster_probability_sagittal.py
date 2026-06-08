#!/usr/bin/env python3
"""Create per-slice MERFISH probability scatter plots for Tangram cluster assignments."""

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sp-h5ad", required=True, help="Spatial AnnData (.h5ad) file")
    parser.add_argument("--metadata-csv", required=True, help="Cell metadata CSV with brain_section_label")
    parser.add_argument("--ccf-csv", required=True, help="CCF coordinates CSV with x,y,z columns")
    parser.add_argument("--cluster-csv", required=True,
                        help="Tangram prediction CSV containing cluster probabilities per cell")
    parser.add_argument("--output-dir", required=True, help="Directory to store slice plots")
    parser.add_argument("--exclude-columns", default="brain_section_label,uniform_density,rna_count_based_density,cluster_id",
                        help="Comma-separated columns to exclude from probability set")
    parser.add_argument("--point-size", type=float, default=2.0, help="Scatter plot point size [default %(default)s]")
    parser.add_argument("--padding", type=float, default=0.05, help="Axis padding in mm [default %(default)s]")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata_sp = sc.read_h5ad(args.sp_h5ad)
    meta = pd.read_csv(args.metadata_csv, index_col="cell_label")
    ccf = pd.read_csv(args.ccf_csv, index_col="cell_label")
    cluster_df = pd.read_csv(args.cluster_csv, index_col=0)

    adata_sp.obs_names = adata_sp.obs_names.astype(str)
    cluster_df.index = cluster_df.index.astype(str)

    shared_cells = adata_sp.obs_names.intersection(meta.index).intersection(ccf.index)
    if not len(shared_cells):
        raise ValueError("No overlapping cells among AnnData, metadata, and CCF tables")

    adata_sp = adata_sp[shared_cells].copy()
    meta = meta.loc[shared_cells]
    ccf = ccf.loc[shared_cells]
    cluster_df = cluster_df.reindex(adata_sp.obs_names).fillna(0)

    adata_sp.obs["x_um"] = ccf["x"]
    adata_sp.obs["y_um"] = ccf["y"]
    adata_sp.obs["z"] = ccf["z"]
    adata_sp.obs["brain_section_label"] = meta["brain_section_label"]

    excluded = set(col.strip() for col in args.exclude_columns.split(",") if col)
    prob_cols = [col for col in cluster_df.columns if col not in excluded]
    if not prob_cols:
        raise ValueError("No probability columns available after exclusions")

    global_min = cluster_df[prob_cols].min().min()
    global_max = cluster_df[prob_cols].max().max()
    x_all = adata_sp.obs["x_um"]
    y_all = adata_sp.obs["y_um"]
    padding = args.padding
    xmin, xmax = x_all.min() - padding, x_all.max() + padding
    ymin, ymax = y_all.min() - padding, y_all.max() + padding

    for cluster in prob_cols:
        cluster_outdir = outdir / f"cluster_{cluster}"
        cluster_outdir.mkdir(parents=True, exist_ok=True)
        adata_sp.obs["prob"] = cluster_df[cluster].reindex(adata_sp.obs_names).fillna(0)

        for slice_name in adata_sp.obs["brain_section_label"].unique():
            df = adata_sp.obs[adata_sp.obs["brain_section_label"] == slice_name].copy()
            if df.empty:
                continue
            df = df.sort_values(by="prob", ascending=True)
            colors = np.full(len(df), "#1B1B1B", dtype=object)
            positive = df["prob"] > 0
            if positive.any():
                scaled = (df.loc[positive, "prob"] - global_min) / (global_max - global_min or 1)
                viridis_colors = plt.cm.viridis(np.clip(scaled, 0, 1))[:, :3]
                colors[positive.to_numpy()] = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in viridis_colors]

            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(df["x_um"], df["y_um"], c=colors, s=args.point_size, alpha=1)
            ax.set_title(f"Cluster {cluster} – {slice_name} (z={df['z'].mean():.2f})")
            ax.set_xlabel("x (Rostral-Caudal, mm)")
            ax.set_ylabel("y (Inferior-Superior, mm)")
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect("equal")
            fig.tight_layout()

            sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=global_min, vmax=global_max))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, orientation="vertical")
            cbar.set_label("Cluster Probability", rotation=270, labelpad=15)

            fig.savefig(cluster_outdir / f"{slice_name}_cluster_{cluster}.pdf")
            plt.close(fig)


if __name__ == "__main__":
    main()
