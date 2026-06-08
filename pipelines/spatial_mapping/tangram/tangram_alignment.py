#!/usr/bin/env python3
"""Run Tangram alignment between single-cell and spatial datasets with configurable options."""

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sc-h5ad", required=True, help="Single-cell reference AnnData (.h5ad)")
    parser.add_argument("--sp-h5ad", required=True, help="Spatial AnnData (.h5ad)")
    parser.add_argument("--gene-map-csv", required=True,
                        help="CSV mapping spatial gene identifiers to symbols")
    parser.add_argument("--gene-id-column", default="gene_identifier",
                        help="Column in gene map representing spatial gene ids [default %(default)s]")
    parser.add_argument("--gene-symbol-column", default="gene_symbol",
                        help="Column in gene map representing gene symbols [default %(default)s]")
    parser.add_argument("--output-dir", required=True, help="Directory for alignment outputs")
    parser.add_argument("--max-cells", type=int, default=10000,
                        help="Maximum number of single-cell observations to map [default %(default)s]")
    parser.add_argument("--hv-genes", type=int, default=750,
                        help="Number of highly variable genes to retain [default %(default)s]")
    parser.add_argument("--epochs", type=int, default=400,
                        help="Number of training epochs [default %(default)s]")
    parser.add_argument("--lambda-reg", type=float, default=0.1,
                        help="Tangram lambda_f_reg parameter [default %(default)s]")
    parser.add_argument("--density-prior", default="uniform",
                        help="Tangram density prior [default %(default)s]")
    parser.add_argument("--mode", default="cells",
                        help="Tangram mode (cells or clusters) [default %(default)s]")
    parser.add_argument("--prefix", default=None,
                        help="Prefix for output files (default derived from spatial filename)")
    parser.add_argument("--no-plots", action="store_true", help="Disable diagnostic plots")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for subsampling")
    return parser.parse_args()


def ensure_output_dir(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)


def load_gene_map(path: Path, id_column: str, symbol_column: str) -> dict:
    mapping_df = pd.read_csv(path)
    if id_column not in mapping_df.columns or symbol_column not in mapping_df.columns:
        raise ValueError(f"Gene map must include '{id_column}' and '{symbol_column}' columns")
    return dict(zip(mapping_df[id_column], mapping_df[symbol_column]))


def plot_residual_variance(adata_sc_tmp, output_dir: Path, prefix: str) -> None:
    mean_expr = adata_sc_tmp.var["means"].values
    residual_var = adata_sc_tmp.var["dispersions_norm"].values

    sorted_idx = np.argsort(residual_var)[::-1]
    num_genes = len(residual_var)
    block_colors = [f"C{i % 10}" for i in range((num_genes // 100) + 1)]
    gene_colors = np.array([block_colors[i // 100] for i in range(num_genes)])[np.argsort(sorted_idx)]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(mean_expr, residual_var, alpha=0.6, s=10, c=gene_colors)
    ax.set_xlabel("Mean Expression (log scale)")
    ax.set_ylabel("Normalised Dispersion")
    ax.set_title("Normalised Dispersion vs Mean Expression")
    ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(output_dir / f"{prefix}_residual_variance_vs_expression.pdf")
    plt.close(fig)


def plot_expression_distributions(adata_sc, adata_sp, output_dir: Path, prefix: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    sc.pl.highest_expr_genes(adata_sc, n_top=20, ax=axes[0], show=False)
    axes[0].set_title("Top 20 Expressed Genes - Single-cell")
    sc.pl.highest_expr_genes(adata_sp, n_top=20, ax=axes[1], show=False)
    axes[1].set_title("Top 20 Expressed Genes - Spatial")
    fig.tight_layout()
    fig.savefig(output_dir / f"{prefix}_expression_distributions.pdf")
    plt.close(fig)


def plot_training_history(history: dict, output_dir: Path, prefix: str, spatial_name: str) -> None:
    for key in ("total_loss", "main_loss", "kl_reg"):
        values = history.get(key)
        if not values:
            continue
        fig, ax = plt.subplots()
        ax.plot([float(v) for v in values])
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Loss")
        ax.set_title(f"Tangram Training: {key}")
        fig.tight_layout()
        fig.savefig(output_dir / f"{spatial_name}_{key}.pdf")
        plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    ensure_output_dir(output_dir)
    prefix = args.prefix or Path(args.sp_h5ad).stem

    np.random.seed(args.seed)

    adata_sc = sc.read_h5ad(args.sc_h5ad)
    adata_sp = sc.read_h5ad(args.sp_h5ad)

    gene_map = load_gene_map(Path(args.gene_map_csv), args.gene_id_column, args.gene_symbol_column)
    adata_sp.var_names = adata_sp.var_names.map(gene_map).astype(str)
    if adata_sp.var_names.str.contains("ENSMUSG").any():
        raise ValueError("Gene name remapping failed; spatial dataset still contains ENSMUSG identifiers")

    shared_genes = sorted(set(adata_sc.var_names).intersection(adata_sp.var_names))
    if not shared_genes:
        raise ValueError("No overlapping genes between single-cell and spatial datasets")

    if "cluster_id" not in adata_sc.obs.columns:
        raise ValueError("Single-cell AnnData must contain a 'cluster_id' column in obs")

    adata_sc = adata_sc[:, shared_genes]
    adata_sp = adata_sp[:, shared_genes]

    if args.max_cells and adata_sc.n_obs > args.max_cells:
        sampled = np.random.choice(adata_sc.obs_names, size=args.max_cells, replace=False)
        adata_sc = adata_sc[sampled].copy()

    adata_sc_tmp = adata_sc[:, shared_genes].copy()
    sc.pp.highly_variable_genes(adata_sc_tmp, flavor="cell_ranger", n_top_genes=None, inplace=True)

    if not args.no_plots:
        plot_residual_variance(adata_sc_tmp, output_dir, prefix)

    sc.pp.highly_variable_genes(adata_sc_tmp, n_top_genes=args.hv_genes, subset=True, flavor="seurat", inplace=True)
    hv_genes = list(adata_sc_tmp.var_names.intersection(shared_genes))
    if not hv_genes:
        raise ValueError("No highly variable genes retained; adjust --hv-genes or preprocessing")
    adata_sc = adata_sc[:, hv_genes]
    adata_sp = adata_sp[:, hv_genes]

    if not args.no_plots:
        plot_expression_distributions(adata_sc, adata_sp, output_dir, prefix)

    tg.pp_adatas(adata_sc, adata_sp, genes=hv_genes)

    mapping_matrix = tg.map_cells_to_space(
        adata_sc=adata_sc,
        adata_sp=adata_sp,
        mode=args.mode,
        density_prior=args.density_prior,
        num_epochs=args.epochs,
        lambda_f_reg=args.lambda_reg
    )

    if not args.no_plots:
        plot_training_history(mapping_matrix.uns.get("training_history", {}), output_dir, prefix, Path(args.sp_h5ad).stem)

    tg.project_cell_annotations(mapping_matrix, adata_sp, annotation="cluster_id")

    if "cluster_id" not in adata_sp.obs.columns:
        if "tangram_ct_pred" not in adata_sp.obsm:
            raise ValueError("Tangram did not produce cluster predictions")
        cluster_probs = adata_sp.obsm["tangram_ct_pred"].to_numpy()
        cluster_ids = cluster_probs.argmax(axis=1)
        labels = adata_sc.obs["cluster_id"].astype(str).unique().tolist()
        adata_sp.obs["cluster_id"] = [labels[i] for i in cluster_ids]

    output_h5ad = output_dir / f"{prefix}_predicted_cluster_id.h5ad"
    output_csv = output_dir / f"{prefix}_predicted_cluster_id.csv"
    adata_sp.write_h5ad(output_h5ad)

    if "tangram_ct_pred" not in adata_sp.obsm:
        adata_sp.obs.to_csv(output_csv)
    else:
        cluster_probs = adata_sp.obsm["tangram_ct_pred"].to_numpy()
        clusters = adata_sc.obs["cluster_id"].astype(str).unique().tolist()
        cluster_probs_df = pd.DataFrame(cluster_probs, index=adata_sp.obs_names, columns=clusters)
        cluster_probs_df.index.name = "cell_id"
        full_metadata = pd.concat([adata_sp.obs, cluster_probs_df], axis=1)
        full_metadata.to_csv(output_csv)

    print(f"Saved: {output_h5ad}")
    print(f"Saved: {output_csv}")


if __name__ == "__main__":
    main()
