"""
run_nichecompass_integrated.py

NicheCompass niche identification for chronic vs healthy wound Xenium data.

Loads converted h5ad files from annotated_data/, integrates all 6 samples,
trains NicheCompass, identifies niches via Leiden clustering, and saves results.

Usage:
    python run_nichecompass_integrated.py [--base_dir ./annotated_data] [--output_dir ./nichecompass_results]

Requires:
    pip install nichecompass scanpy squidpy anndata
    (with PyTorch + PyG already installed for GPU support)
"""

import os
import sys
import glob
import argparse
import warnings
from datetime import datetime

# Force non-interactive backend BEFORE importing pyplot anywhere
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import scipy.sparse as sp

# Suppress scanpy/squidpy interactive plot attempts
sc.settings.autoshow = False

from nichecompass.models import NicheCompass
from nichecompass.utils import (
    add_gps_from_gp_dict_to_adata,
    create_new_color_dict,
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
    filter_and_combine_gp_dict_gps_v2,
)

warnings.filterwarnings("ignore")


# =============================================================================
# CONFIG — edit these if needed
# =============================================================================

# Data keys (matching your h5ad structure)
SPATIAL_KEY = "X_spatial"          # obsm key for proseg centroid coordinates
CELL_TYPE_KEY = "predicted_cell_type"  # obs column for cell annotations
COUNTS_KEY = "counts"              # layer key NicheCompass expects for raw counts

# Condition mapping — prefix -> condition label
CONDITION_MAP = {
    "HS": "healthy",
    "LS": "chronic",
}

# Spatial graph
N_NEIGHBORS = 8  # User guide: smaller = better NID; 8 is good for single-cell

# Gene programs
SPECIES = "human"

# NicheCompass AnnData keys
ADJ_KEY = "spatial_connectivities"
GP_NAMES_KEY = "nichecompass_gp_names"
ACTIVE_GP_NAMES_KEY = "nichecompass_active_gp_names"
GP_TARGETS_MASK_KEY = "nichecompass_gp_targets"
GP_TARGETS_CATEGORIES_MASK_KEY = "nichecompass_gp_targets_categories"
GP_SOURCES_MASK_KEY = "nichecompass_gp_sources"
GP_SOURCES_CATEGORIES_MASK_KEY = "nichecompass_gp_sources_categories"
LATENT_KEY = "nichecompass_latent"

# Architecture — GATv2 recommended for single-cell resolution data
CONV_LAYER_ENCODER = "gatv2conv"
ACTIVE_GP_THRESH_RATIO = 0.01

# Training
N_EPOCHS = 400
N_EPOCHS_ALL_GPS = 25
LR = 0.0001
LAMBDA_EDGE_RECON = 50000.0       
LAMBDA_GENE_EXPR_RECON = 300.0
LAMBDA_L1_MASKED = 0.0
LAMBDA_L1_ADDON = 30.0
EDGE_BATCH_SIZE = 512             
N_SAMPLED_NEIGHBORS = 4

# Clustering
LATENT_LEIDEN_RESOLUTION = 0.8
LATENT_CLUSTER_KEY = "nichecompass_niche"


# =============================================================================
# STEP 1: Load and prepare data
# =============================================================================

def discover_h5ad_files(base_dir):
    """Find all _annotated.h5ad files and assign sample IDs and conditions."""
    pattern = os.path.join(base_dir, "*", "*_annotated.h5ad")
    files = sorted(glob.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No *_annotated.h5ad files found under {base_dir}")

    samples = []
    for f in files:
        # Extract sample name from parent directory (e.g., "HS_Region_2")
        sample_id = os.path.basename(os.path.dirname(f))
        prefix = sample_id.split("_")[0]  # "HS" or "LS"
        condition = CONDITION_MAP.get(prefix)
        if condition is None:
            print(f"[WARN] Unknown prefix '{prefix}' for {sample_id}, skipping")
            continue
        samples.append({
            "sample_id": sample_id,
            "condition": condition,
            "path": f,
        })

    print(f"Discovered {len(samples)} samples:")
    for s in samples:
        print(f"  {s['sample_id']} ({s['condition']}): {s['path']}")

    return samples


def load_and_prepare_sample(sample_info):
    """
    Load a single h5ad, add metadata, prepare counts layer and spatial graph.
    """
    sample_id = sample_info["sample_id"]
    condition = sample_info["condition"]
    path = sample_info["path"]

    print(f"\n[LOAD] {sample_id} from {path}")
    adata = sc.read_h5ad(path)

    # --- Add sample/condition metadata ---
    adata.obs["sample_id"] = sample_id
    adata.obs["condition"] = condition

    # --- Prepare counts layer ---
    # scCustomize put raw counts into .X; NicheCompass expects a named layer.
    # CRITICAL: cast to float32 — NicheCompass / PyTorch expects float32,
    # but h5ad files often store .X as float64 which causes dtype mismatch
    # errors in the encoder ("mat1 and mat2 must have the same dtype").
    if COUNTS_KEY not in adata.layers:
        if sp.issparse(adata.X):
            adata.layers[COUNTS_KEY] = adata.X.astype(np.float32).copy()
        else:
            adata.layers[COUNTS_KEY] = adata.X.astype(np.float32).copy()
        print(f"  Copied .X to .layers['{COUNTS_KEY}'] (cast to float32)")
    else:
        # Ensure existing counts layer is also float32
        if sp.issparse(adata.layers[COUNTS_KEY]):
            if adata.layers[COUNTS_KEY].dtype != np.float32:
                adata.layers[COUNTS_KEY] = adata.layers[COUNTS_KEY].astype(np.float32)
                print(f"  Cast existing .layers['{COUNTS_KEY}'] to float32")
        else:
            if adata.layers[COUNTS_KEY].dtype != np.float32:
                adata.layers[COUNTS_KEY] = adata.layers[COUNTS_KEY].astype(np.float32)
                print(f"  Cast existing .layers['{COUNTS_KEY}'] to float32")

    # Also ensure .X itself is float32 for consistency
    if sp.issparse(adata.X):
        if adata.X.dtype != np.float32:
            adata.X = adata.X.astype(np.float32)
    else:
        if adata.X.dtype != np.float32:
            adata.X = adata.X.astype(np.float32)

    # --- Ensure spatial coordinates exist ---
    if SPATIAL_KEY not in adata.obsm:
        raise KeyError(f"Spatial key '{SPATIAL_KEY}' not in obsm for {sample_id}")

    # --- Build per-sample spatial graph ---
    # Must be done BEFORE concatenation so samples remain disconnected
    sq.gr.spatial_neighbors(
        adata,
        spatial_key=SPATIAL_KEY,
        n_neighs=N_NEIGHBORS,
        coord_type="generic",
    )
    n_edges = adata.obsp["spatial_connectivities"].nnz
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes, {n_edges} spatial edges")

    # --- Ensure var_names are unique ---
    adata.var_names_make_unique()

    # --- Make obs index unique by prefixing sample_id ---
    adata.obs_names = [f"{sample_id}_{i}" for i in adata.obs_names]

    return adata


def concatenate_samples(adatas):
    """
    Concatenate per-sample AnnDatas with block-diagonal spatial graph.
    """
    from scipy.sparse import block_diag

    print(f"\n[CONCAT] Merging {len(adatas)} samples...")

    # Collect spatial graphs before concat
    connectivities = [a.obsp["spatial_connectivities"] for a in adatas]
    distances = [a.obsp["spatial_distances"] for a in adatas]

    # Concatenate — use inner join on var to handle slight gene differences
    adata = ad.concat(adatas, join="inner", merge="same")

    # Rebuild block-diagonal spatial graph (samples stay disconnected)
    adata.obsp["spatial_connectivities"] = block_diag(connectivities, format="csr")
    adata.obsp["spatial_distances"] = block_diag(distances, format="csr")

    # Ensure categorical types
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
    adata.obs["condition"] = adata.obs["condition"].astype("category")
    adata.obs[CELL_TYPE_KEY] = adata.obs[CELL_TYPE_KEY].astype("category")

    print(f"  Combined: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"  Samples: {adata.obs['sample_id'].value_counts().to_dict()}")
    print(f"  Conditions: {adata.obs['condition'].value_counts().to_dict()}")

    return adata


# =============================================================================
# STEP 2: Build gene program prior mask
# =============================================================================

def build_gene_programs(adata):
    """
    Extract and combine gene program (GP) databases for the prior mask.
    NicheCompass uses these to make its latent space interpretable.

    MEBOCOST requires bundled data files that may not be present outside
    the NicheCompass repo. If it fails, we proceed with OmniPath + NicheNet
    which already provide strong ligand-receptor coverage.
    """
    print("\n[GP] Building gene program prior mask...")

    # Extract from three databases
    print("  Extracting OmniPath LR interactions...")
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=SPECIES,
        min_curation_effort=0,
    )

    print("  Extracting NicheNet LRT interactions...")
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=SPECIES,
        version="v2",
        keep_target_genes_ratio=1.,
        max_n_target_genes_per_gp=250,
    )

    mebocost_gp_dict = None
    try:
        print("  Extracting MEBOCOST MS interactions...")
        mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
            species=SPECIES,
        )
    except (FileNotFoundError, Exception) as e:
        print(f"  [WARN] MEBOCOST extraction failed: {e}")
        print("  Proceeding with OmniPath + NicheNet only.")

    # Combine and filter — function takes a positional list of GP dicts,
    # NOT named keyword arguments
    print("  Combining and filtering gene programs...")
    gp_dicts = [omnipath_gp_dict, nichenet_gp_dict]
    if mebocost_gp_dict is not None:
        gp_dicts.append(mebocost_gp_dict)

    combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
        gp_dicts,
        verbose=True,
    )

    # Add GPs to adata — note: active_gp_names_key is NOT a parameter of
    # this function (it is only used at NicheCompass model init)
    add_gps_from_gp_dict_to_adata(
        gp_dict=combined_gp_dict,
        adata=adata,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        gp_names_key=GP_NAMES_KEY,
        min_genes_per_gp=2,
        min_source_genes_per_gp=1,
        min_target_genes_per_gp=1,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None,
    )

    n_gps = len(adata.uns.get(GP_NAMES_KEY, []))
    print(f"  {n_gps} total gene programs")

    return adata


# =============================================================================
# STEP 3: Train NicheCompass
# =============================================================================

def train_nichecompass(adata):
    """Initialize and train the NicheCompass model."""
    print("\n[TRAIN] Initializing NicheCompass model...")

    model = NicheCompass(
        adata,
        counts_key=COUNTS_KEY,
        adj_key=ADJ_KEY,
        cat_covariates_embeds_injection=["gene_expr_decoder"],
        cat_covariates_keys=["sample_id"],
        cat_covariates_no_edges=[True],
        cat_covariates_embeds_nums=[len(adata.obs["sample_id"].cat.categories)],
        gp_names_key=GP_NAMES_KEY,
        active_gp_names_key=ACTIVE_GP_NAMES_KEY,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        latent_key=LATENT_KEY,
        conv_layer_encoder=CONV_LAYER_ENCODER,
        active_gp_thresh_ratio=ACTIVE_GP_THRESH_RATIO,
    )

    print("[TRAIN] Starting training...")
    model.train(
        n_epochs=N_EPOCHS,
        n_epochs_all_gps=N_EPOCHS_ALL_GPS,
        lr=LR,
        lambda_edge_recon=LAMBDA_EDGE_RECON,
        lambda_gene_expr_recon=LAMBDA_GENE_EXPR_RECON,
        lambda_l1_masked=LAMBDA_L1_MASKED,
        lambda_l1_addon=LAMBDA_L1_ADDON,
        edge_batch_size=EDGE_BATCH_SIZE,
        n_sampled_neighbors=N_SAMPLED_NEIGHBORS,
        use_cuda_if_available=True,
        verbose=False,
    )

    print("[TRAIN] Training complete.")
    return model


# =============================================================================
# STEP 4: Cluster niches and compute embeddings
# =============================================================================

def cluster_niches(model, resolution=LATENT_LEIDEN_RESOLUTION):
    """
    Compute neighbor graph on latent space, run Leiden clustering, compute UMAP.
    """
    print(f"\n[CLUSTER] Leiden clustering at resolution={resolution}...")

    adata = model.adata

    # Neighbor graph on NicheCompass latent space
    sc.pp.neighbors(
        adata,
        use_rep=LATENT_KEY,
        key_added=LATENT_KEY,
    )

    # Leiden clustering
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=LATENT_CLUSTER_KEY,
        neighbors_key=LATENT_KEY,
    )

    # UMAP for visualization
    sc.tl.umap(adata, neighbors_key=LATENT_KEY)

    n_niches = adata.obs[LATENT_CLUSTER_KEY].nunique()
    print(f"  Identified {n_niches} niches")

    # Summarize niche composition
    ct_by_niche = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs[CELL_TYPE_KEY],
        normalize="index",
    )
    print("\n  Top cell types per niche:")
    for niche in ct_by_niche.index:
        top3 = ct_by_niche.loc[niche].nlargest(3)
        top_str = ", ".join([f"{ct} ({pct:.0%})" for ct, pct in top3.items()])
        print(f"    Niche {niche}: {top_str}")

    return adata


# =============================================================================
# STEP 5: Save results
# =============================================================================

def save_results(adata, model, output_dir):
    """
    Save all outputs in organized subdirectories:
        output_dir/
            objects/         — h5ad and model
            tables/          — CSV summary tables
            plots/
                umap/        — UMAP embeddings
                spatial/     — per-sample spatial niche maps
                composition/ — niche composition heatmaps/bars
    """
    import matplotlib.pyplot as plt

    # --- Create directory structure ---
    dirs = {
        "objects":     os.path.join(output_dir, "objects"),
        "tables":      os.path.join(output_dir, "tables"),
        "plots_umap":  os.path.join(output_dir, "plots", "umap"),
        "plots_spatial": os.path.join(output_dir, "plots", "spatial"),
        "plots_composition": os.path.join(output_dir, "plots", "composition"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    samples = adata.obs["sample_id"].cat.categories.tolist()

    # ---- OBJECTS ----
    adata_path = os.path.join(dirs["objects"], "nichecompass_integrated.h5ad")
    print(f"\n[SAVE] Writing {adata_path}")
    adata.write_h5ad(adata_path)

    model_dir = os.path.join(dirs["objects"], f"model_{timestamp}")
    print(f"[SAVE] Saving model to {model_dir}")
    model.save(dir_path=model_dir, overwrite=True, save_adata=False)

    # ---- TABLES ----
    # Per-cell niche assignments
    niche_df = adata.obs[
        ["sample_id", "condition", CELL_TYPE_KEY, LATENT_CLUSTER_KEY]
    ].copy()
    niche_csv = os.path.join(dirs["tables"], "niche_assignments.csv")
    niche_df.to_csv(niche_csv)
    print(f"[SAVE] {niche_csv}")

    # Niche x cell type composition (row-normalized)
    ct_by_niche = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs[CELL_TYPE_KEY],
        normalize="index",
    )
    ct_by_niche.to_csv(os.path.join(dirs["tables"], "niche_composition.csv"))

    # Niche x sample frequency (column-normalized)
    niche_by_sample = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs["sample_id"],
        normalize="columns",
    )
    niche_by_sample.to_csv(os.path.join(dirs["tables"], "niche_sample_frequency.csv"))

    # Niche x condition frequency (column-normalized)
    niche_by_cond = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs["condition"],
        normalize="columns",
    )
    niche_by_cond.to_csv(os.path.join(dirs["tables"], "niche_condition_frequency.csv"))

    # Niche cell counts (absolute)
    niche_counts = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs["sample_id"],
    )
    niche_counts.to_csv(os.path.join(dirs["tables"], "niche_cell_counts.csv"))

    print(f"[SAVE] All tables -> {dirs['tables']}")

    # ---- PLOTS: UMAP ----
    print("[PLOT] Generating UMAP plots...")

    # UMAP colored by niche
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=LATENT_CLUSTER_KEY, ax=ax, show=False,
               title="NicheCompass Niches", legend_loc="on data")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_niches.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # UMAP colored by sample
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="sample_id", ax=ax, show=False,
               title="Samples in Latent Space")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_samples.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # UMAP colored by condition
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="condition", ax=ax, show=False,
               title="Condition in Latent Space")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_condition.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # UMAP colored by cell type
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.umap(adata, color=CELL_TYPE_KEY, ax=ax, show=False,
               title="Cell Types in Latent Space")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_cell_types.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    print(f"  UMAP plots -> {dirs['plots_umap']}")

    # ---- PLOTS: SPATIAL per sample ----
    print("[PLOT] Generating spatial niche maps per sample...")

    for sample in samples:
        adata_sub = adata[adata.obs["sample_id"] == sample].copy()
        condition = adata_sub.obs["condition"].iloc[0]

        # Spatial scatter colored by niche
        fig, ax = plt.subplots(figsize=(10, 10))
        coords = adata_sub.obsm[SPATIAL_KEY]
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=adata_sub.obs[LATENT_CLUSTER_KEY].cat.codes,
            cmap="tab20", s=1, alpha=0.8, rasterized=True,
        )
        ax.set_title(f"{sample} ({condition}) — Niches")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        plt.colorbar(scatter, ax=ax, label="Niche")
        fig.savefig(
            os.path.join(dirs["plots_spatial"], f"{sample}_niches.pdf"),
            bbox_inches="tight", dpi=150,
        )
        plt.close(fig)

        # Spatial scatter colored by cell type
        fig, ax = plt.subplots(figsize=(10, 10))
        ct_codes = adata_sub.obs[CELL_TYPE_KEY].cat.codes
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=ct_codes, cmap="tab20", s=1, alpha=0.8, rasterized=True,
        )
        ax.set_title(f"{sample} ({condition}) — Cell Types")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        fig.savefig(
            os.path.join(dirs["plots_spatial"], f"{sample}_cell_types.pdf"),
            bbox_inches="tight", dpi=150,
        )
        plt.close(fig)

    print(f"  Spatial plots -> {dirs['plots_spatial']}")

    # ---- PLOTS: COMPOSITION ----
    print("[PLOT] Generating composition plots...")

    # Heatmap: niche x cell type composition
    fig, ax = plt.subplots(figsize=(max(12, ct_by_niche.shape[1] * 0.5),
                                    max(6, ct_by_niche.shape[0] * 0.4)))
    import seaborn as sns
    sns.heatmap(ct_by_niche, cmap="YlOrRd", ax=ax, linewidths=0.5,
                xticklabels=True, yticklabels=True)
    ax.set_title("Cell Type Composition per Niche")
    ax.set_ylabel("Niche")
    ax.set_xlabel("Cell Type")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    fig.savefig(
        os.path.join(dirs["plots_composition"], "niche_celltype_heatmap.pdf"),
        bbox_inches="tight", dpi=150,
    )
    plt.close(fig)

    # Stacked bar: niche proportions per sample, grouped by condition
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for idx, cond in enumerate(["healthy", "chronic"]):
        cond_samples = [s for s in samples if adata.obs.loc[
            adata.obs["sample_id"] == s, "condition"].iloc[0] == cond]
        if not cond_samples:
            continue
        cond_freq = niche_by_sample[cond_samples]
        cond_freq.T.plot(kind="bar", stacked=True, ax=axes[idx],
                         colormap="tab20", legend=False)
        axes[idx].set_title(f"{cond.capitalize()} Wounds")
        axes[idx].set_xlabel("Sample")
        axes[idx].set_ylabel("Niche Proportion")
        axes[idx].tick_params(axis="x", rotation=45)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Niche", bbox_to_anchor=(1.02, 0.5),
               loc="center left", fontsize=8)
    fig.suptitle("Niche Proportions per Sample by Condition", y=1.02)
    fig.savefig(
        os.path.join(dirs["plots_composition"], "niche_proportions_by_condition.pdf"),
        bbox_inches="tight", dpi=150,
    )
    plt.close(fig)

    # Bar: niche proportions averaged per condition (side-by-side)
    fig, ax = plt.subplots(figsize=(max(8, niche_by_cond.shape[0] * 0.6), 5))
    niche_by_cond.plot(kind="bar", ax=ax)
    ax.set_title("Niche Proportions by Condition")
    ax.set_xlabel("Niche")
    ax.set_ylabel("Proportion of Cells")
    ax.legend(title="Condition")
    plt.xticks(rotation=0)
    fig.savefig(
        os.path.join(dirs["plots_composition"], "niche_proportions_condition_comparison.pdf"),
        bbox_inches="tight", dpi=150,
    )
    plt.close(fig)

    print(f"  Composition plots -> {dirs['plots_composition']}")
    print(f"\n[DONE] All outputs in {output_dir}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="NicheCompass niche analysis")
    parser.add_argument(
        "--base_dir", default="./annotated_data",
        help="Directory containing sample subfolders with .h5ad files",
    )
    parser.add_argument(
        "--output_dir", default="./nichecompass_results",
        help="Output directory for results",
    )
    parser.add_argument(
        "--leiden_resolution", type=float, default=LATENT_LEIDEN_RESOLUTION,
        help="Leiden clustering resolution (default: 0.8)",
    )
    args = parser.parse_args()

    print("=" * 70)
    print("NicheCompass Pipeline — Chronic vs Healthy Wound Niche Analysis")
    print("=" * 70)

    # Step 1: Discover and load samples
    sample_info_list = discover_h5ad_files(args.base_dir)
    adatas = [load_and_prepare_sample(s) for s in sample_info_list]
    adata = concatenate_samples(adatas)

    # Free per-sample objects
    del adatas

    # Step 2: Build gene programs
    adata = build_gene_programs(adata)

    # Step 3: Train model
    model = train_nichecompass(adata)

    # Step 4: Cluster niches
    adata = cluster_niches(model, resolution=args.leiden_resolution)

    # Step 5: Save results
    save_results(adata, model, args.output_dir)


if __name__ == "__main__":
    main()