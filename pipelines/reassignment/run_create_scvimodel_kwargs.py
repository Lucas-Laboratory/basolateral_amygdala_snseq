#!/usr/bin/env python3
"""Train or reuse SCVI/SCANVI models for barcode reassignment workflows with configurable I/O."""
from __future__ import annotations

import argparse
import os
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from scipy.stats import zscore

try:
    from torchviz import make_dot
except ImportError:  # pragma: no cover - optional dependency
    make_dot = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-h5ad", required=True, help="Training AnnData (.h5ad) with cluster labels")
    parser.add_argument("--target-h5ad", required=True, help="Target AnnData (.h5ad) for prediction")
    parser.add_argument("--model-dir", default=".", help="Directory to store/load SCVI/SCANVI models")
    parser.add_argument("--output-dir", default="scanvi_outputs", help="Directory for plots and activation CSVs")
    parser.add_argument("--predictions-csv", default=None,
                        help="Path for predicted clusters CSV [default: <output-dir>/scanvi_target_predictions.csv]")
    parser.add_argument("--retrain-scvi", action="store_true", help="Force retraining of SCVI")
    parser.add_argument("--retrain-scanvi", action="store_true", help="Force retraining of SCANVI")
    parser.add_argument("--predict-only", action="store_true",
                        help="Skip training/prediction steps (useful to just validate inputs)")
    parser.add_argument("--skip-confirm", action="store_true",
                        help="Skip interactive confirmation before clipping HVGs")
    parser.add_argument("--export-csvs", action="store_true",
                        help="Export activation CSVs per cluster (train/target)")
    parser.add_argument("--suppress-csv", action="store_true",
                        help="Suppress target activation CSV export even if --export-csvs set")
    parser.add_argument("--num-threads", type=int, default=8, help="Number of BLAS threads to use [default %(default)s]")
    parser.add_argument("--max-epochs-scvi", type=int, default=300, help="Maximum epochs for SCVI training")
    parser.add_argument("--max-epochs-scanvi", type=int, default=150, help="Maximum epochs for SCANVI training")
    parser.add_argument("--latent-dim", type=int, default=30, help="Latent dimensionality for SCVI")
    return parser.parse_args()


def configure_threads(num_threads: int) -> None:
    for env in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[env] = str(num_threads)


def ensure_predictions_header(path: Path) -> None:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("barcode,predicted_cluster\n", encoding="utf-8")


def report_matrix_health(name: str, matrix) -> None:
    print(f"=== Checking for NaN and Inf in {name}.X ===")
    if isinstance(matrix, np.ndarray):
        print("Dense matrix detected.")
        print("NaN count:", np.isnan(matrix).sum())
        print("Inf count:", np.isinf(matrix).sum())
    elif hasattr(matrix, "data"):
        print("Sparse matrix detected.")
        print("NaN count:", np.isnan(matrix.data).sum())
        print("Inf count:", np.isinf(matrix.data).sum())
    else:
        print("Unknown matrix type.")

    print(f"=== Inspecting extreme values in {name}.X ===")
    if hasattr(matrix, "data"):
        print("Max value:", matrix.data.max())
        print("Min value:", matrix.data.min())
    else:
        print("Max value:", matrix.max())
        print("Min value:", matrix.min())


def clip_matrix(matrix, upper: float = 1e6):
    if hasattr(matrix, "data"):
        matrix.data = np.clip(matrix.data, 0, upper)
    else:
        matrix[:] = np.clip(matrix, 0, upper)


def save_curve(fig_path: Path, history_key: str, history: dict, title: str, ylabel: str, legend_label: str) -> None:
    values = history.get(history_key)
    if not values:
        return
    plt.figure()
    plt.plot(values, label=legend_label)
    val_key = f"{history_key}_validation"
    if history.get(val_key) is not None:
        plt.plot(history[val_key], label=f"Validation {ylabel}")
    plt.xlabel("Epoch")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_path)
    plt.close()


def maybe_render_network(vae: scvi.model.SCVI, output_dir: Path) -> None:
    if make_dot is None:
        print("torchviz not installed; skipping network diagram")
        return
    try:
        vae.module.eval()
        dummy_input = torch.randn((1, vae.module.input_dim)).to(vae.device)
        encoder_out = vae.module.encoder(dummy_input)
        dot = make_dot(encoder_out, params=dict(vae.module.encoder.named_parameters()))
        dot.render(str(output_dir / "scvi_network"), format="pdf")
    except Exception as exc:  # pragma: no cover - best effort
        print(f"Failed to generate SCVI network diagram: {exc}")


def main() -> None:
    args = parse_args()
    configure_threads(args.num_threads)
    warnings.filterwarnings("ignore", category=SyntaxWarning)

    output_dir = Path(args.output_dir)
    activation_dir = output_dir / "scanvi_activation_heatmaps"
    activation_dir.mkdir(parents=True, exist_ok=True)

    models_dir = Path(args.model_dir)
    models_dir.mkdir(parents=True, exist_ok=True)
    scvi_model_path = models_dir / "trained_scvi_model"
    scanvi_model_path = models_dir / "trained_scanvi_model"

    predictions_csv = Path(args.predictions_csv) if args.predictions_csv else output_dir / "scanvi_target_predictions.csv"
    ensure_predictions_header(predictions_csv)

    train_adata = sc.read_h5ad(args.train_h5ad)
    target_adata = sc.read_h5ad(args.target_h5ad)

    report_matrix_health("train_adata", train_adata.X)
    if not args.skip_confirm:
        proceed = input("Proceed with clipping and HVG calculation? (y/n): ")
        if proceed.lower() != "y":
            print("Aborting as requested.")
            return
    if hasattr(train_adata.X, "sort_indices") and not train_adata.X.has_sorted_indices:
        train_adata.X.sort_indices()
    clip_matrix(train_adata.X)

    sc.pp.highly_variable_genes(train_adata, n_top_genes=2000, subset=True, flavor="cell_ranger")

    report_matrix_health("target_adata", target_adata.X)
    clip_matrix(target_adata.X)

    train_adata.obs["labels"] = train_adata.obs["seurat_clusters"].astype(str)
    scvi.model.SCVI.setup_anndata(train_adata, labels_key="labels")

    scvi_exists = scvi_model_path.exists()
    scanvi_exists = scanvi_model_path.exists()

    if args.retrain_scvi or not scvi_exists:
        print("Training SCVI model...")
        vae = scvi.model.SCVI(train_adata, n_latent=args.latent_dim)
        plan_kwargs = {"kl_weight": 1e-4, "n_epochs_kl_warmup": 150, "weight_decay": 1e-4}
        vae.train(max_epochs=args.max_epochs_scvi, train_size=0.8, accelerator="gpu",
                  early_stopping=True, early_stopping_patience=35, plan_kwargs=plan_kwargs)
        vae.save(scvi_model_path, overwrite=True)
        history = vae.history
        save_curve(activation_dir / "scvi_elbo_curve.pdf", "elbo_train", history,
                   "SCVI Training and Validation ELBO", "ELBO", "Train ELBO")
        save_curve(activation_dir / "scvi_reconstruction_loss.pdf", "reconstruction_loss_train", history,
                   "SCVI Reconstruction Loss", "Loss", "Train Reconstruction Loss")
        save_curve(activation_dir / "scvi_kl_divergence.pdf", "kl_local_train", history,
                   "SCVI KL Divergence", "KL Divergence", "Train KL Divergence")
        maybe_render_network(vae, activation_dir)
    else:
        print("Loading pretrained SCVI model...")
        vae = scvi.model.SCVI.load(scvi_model_path, train_adata)

    if args.predict_only:
        print("Predict-only flag supplied; skipping SCANVI training and inference.")
        return

    if args.retrain_scanvi or (not scanvi_exists and scvi_exists):
        print("Training SCANVI model...")
        model = scvi.model.SCANVI.from_scvi_model(vae, labels_key="labels", unlabeled_category="Unknown")
        model.train(max_epochs=args.max_epochs_scanvi, accelerator="gpu",
                    early_stopping=True, early_stopping_patience=20)
        model.save(scanvi_model_path, overwrite=True)
        history = model.history
        save_curve(activation_dir / "scanvi_elbo_curve.pdf", "elbo_train", history,
                   "SCANVI Training and Validation ELBO", "ELBO", "Train ELBO")
        save_curve(activation_dir / "scanvi_reconstruction_loss.pdf", "reconstruction_loss_train", history,
                   "SCANVI Reconstruction Loss", "Loss", "Train Reconstruction Loss")
        save_curve(activation_dir / "scanvi_kl_divergence.pdf", "kl_local_train", history,
                   "SCANVI KL Divergence", "KL Divergence", "Train KL Divergence")
    else:
        print("Loading pretrained SCANVI model...")
        model = scvi.model.SCANVI.load(scanvi_model_path, train_adata)

    print(model.module)

    target_adata.obs["labels"] = "Unknown"
    target_adata = target_adata[:, train_adata.var_names].copy()
    combined = train_adata.concatenate(target_adata, batch_key="dataset", batch_categories=["train", "target"])
    model = scvi.model.SCANVI.load_query_data(combined, model)

    predictions = model.predict(combined)
    prediction_probs = model.predict(combined, soft=True)
    target_mask = combined.obs["dataset"] == "target"
    target_preds = predictions[target_mask]
    target_cells = combined.obs_names[target_mask]
    prob_df = prediction_probs[target_mask].copy()
    prob_df.columns = [f"prob_cluster_{i}" for i in range(prob_df.shape[1])]
    prediction_df = pd.concat([
        pd.DataFrame({"barcode": target_cells, "predicted_cluster": target_preds}).reset_index(drop=True),
        prob_df.reset_index(drop=True)
    ], axis=1)
    prediction_df.to_csv(predictions_csv, index=False)
    print(f"Prediction complete. Output: {predictions_csv}")

    activations: dict[str, list[np.ndarray]] = defaultdict(list)
    barcode_tracker: dict[str, list[str]] = defaultdict(list)

    def register_hooks(scanvi_model: scvi.model.SCANVI) -> None:
        def make_hook(layer_name: str):
            def hook_fn(module, _, output):
                out = output[0] if isinstance(output, tuple) else output
                if isinstance(out, torch.Tensor):
                    activations[layer_name].append(out.detach().cpu().numpy())
                    barcode_tracker[layer_name].extend(current_barcodes)
            return hook_fn

        targets = {
            "latent_z_vector": scanvi_model.module.z_encoder.mean_encoder,
            "linear_128_pre_activation": scanvi_model.module.classifier.classifier[0].fc_layers[0][0],
            "linear_128_post_transform": scanvi_model.module.classifier.classifier[0].fc_layers[0][3],
            "cluster_logits": scanvi_model.module.classifier.classifier[1],
        }
        for name, layer in targets.items():
            layer.register_forward_hook(make_hook(name))

    current_barcodes: list[str] = []
    train_loader = model._make_data_loader(train_adata, batch_size=128)
    register_hooks(model)
    model.module.eval()
    offset = 0
    with torch.no_grad():
        for batch in train_loader:
            batch_x = batch["X"].to(model.device)
            batch_size = batch_x.shape[0]
            batch_ids = list(train_adata.obs_names[offset:offset + batch_size])
            current_barcodes[:] = batch_ids
            offset += batch_size
            model.module.classify(batch_x)

    if args.export_csvs:
        train_dir = activation_dir / "interpretable_train"
        train_dir.mkdir(parents=True, exist_ok=True)
        for layer_name, tensor_list in activations.items():
            try:
                data = np.concatenate(tensor_list, axis=0)
            except Exception as exc:  # pragma: no cover - defensive
                print(f"Skipping {layer_name} (train): {exc}")
                continue
            df = pd.DataFrame(data, index=barcode_tracker[layer_name])
            df["cluster"] = df.index.map(lambda b: train_adata.obs.loc[b, "labels"] if b in train_adata.obs.index else np.nan)
            if df["cluster"].isna().all():
                continue
            z_scores = df.drop(columns=["cluster"]).apply(zscore, axis=0)
            z_scores.insert(0, "cluster", df["cluster"])
            for cl, grp in z_scores.groupby("cluster"):
                mean_row = grp.drop(columns=["cluster"]).mean(axis=0)
                mean_row["cluster"] = f"mean_{cl}"
                out_df = pd.concat([grp, pd.DataFrame([mean_row])], ignore_index=True)
                out_df.to_csv(train_dir / f"train_cluster{cl}_{layer_name}.csv", index=False)

    activations.clear()
    barcode_tracker.clear()
    target_loader = model._make_data_loader(target_adata, batch_size=128)
    register_hooks(model)
    current_barcodes = []
    offset = 0
    with torch.no_grad():
        for batch in target_loader:
            batch_x = batch["X"].to(model.device)
            batch_size = batch_x.shape[0]
            batch_ids = list(target_adata.obs_names[offset:offset + batch_size])
            current_barcodes[:] = batch_ids
            offset += batch_size
            model.module.classify(batch_x)

    if args.export_csvs and not args.suppress_csv:
        target_dir = activation_dir / "interpretable_target"
        target_dir.mkdir(parents=True, exist_ok=True)
        predictions_df = pd.read_csv(predictions_csv)
        predictions_df["barcode"] = predictions_df["barcode"].str.replace("-target$", "", regex=True)
        cluster_lookup = dict(zip(predictions_df["barcode"], predictions_df["predicted_cluster"]))

        for layer_name, tensor_list in activations.items():
            try:
                data = np.concatenate(tensor_list, axis=0)
            except Exception as exc:
                print(f"Skipping {layer_name} (target): {exc}")
                continue
            df = pd.DataFrame(data, index=barcode_tracker[layer_name])
            df.index = df.index.str.replace("-target$", "", regex=True)
            df = df.loc[df.index.intersection(cluster_lookup.keys())]
            if df.empty:
                continue
            df["cluster"] = df.index.map(cluster_lookup.get)
            if df["cluster"].isna().all():
                continue
            z_scores = df.drop(columns=["cluster"]).apply(zscore, axis=0)
            z_scores.insert(0, "cluster", df["cluster"])
            for cl, grp in z_scores.groupby("cluster"):
                mean_row = grp.drop(columns=["cluster"]).mean(axis=0)
                mean_row["cluster"] = f"mean_{cl}"
                out_df = pd.concat([grp, pd.DataFrame([mean_row])], ignore_index=True)
                out_df.to_csv(target_dir / f"target_cluster{cl}_{layer_name}.csv", index=False)

    print("Processing complete.")


if __name__ == "__main__":
    main()
