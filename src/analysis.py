"""
Composition, senescence, UMAP, and volcano figures plus summary tables for the app.

All summaries are small CSVs so the dashboard runs without the large h5ad.
"""

import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pseudobulk import get_counts

# Well-established senescence and SASP markers (mouse symbols).
SENESCENCE_GENES = ["Cdkn2a", "Cdkn1a", "Trp53", "Glb1", "Serpine1",
                    "Il6", "Cxcl1", "Mmp3", "Cdkn2b"]
COLORS = {"young": "#2c7fb8", "old": "#d95f0e"}


def composition(adata, results_dir, data_dir):
    obs = adata.obs
    counts = obs.groupby(["age_group", "_mouse", "_celltype"]).size().rename("n").reset_index()
    totals = counts.groupby(["age_group", "_mouse"])["n"].transform("sum")
    counts["frac"] = counts["n"] / totals
    comp = counts.groupby(["age_group", "_celltype"])["frac"].mean().reset_index()
    comp.to_csv(os.path.join(data_dir, "composition.csv"), index=False)

    pivot = comp.pivot(index="_celltype", columns="age_group", values="frac").fillna(0)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.45 * len(pivot))))
    pivot.plot(kind="barh", ax=ax, color=COLORS)
    ax.set_xlabel("fraction of cells")
    ax.set_ylabel("")
    ax.set_title("Cell-type composition: young vs old")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "composition.png"), dpi=150)
    plt.close(fig)
    return comp


def _lognorm(X):
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    libsize = np.asarray(X.sum(axis=1)).ravel()
    libsize[libsize == 0] = 1.0
    Xn = X.multiply((1e4 / libsize)[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn


def senescence(adata, results_dir, data_dir):
    X, var = get_counts(adata)
    var_index = {g: i for i, g in enumerate(var)}
    present = [g for g in SENESCENCE_GENES if g in var_index]
    if not present:
        print("[senescence] no marker genes present, skipping")
        return None
    Xn = _lognorm(X)
    cols = [var_index[g] for g in present]
    adata.obs["senescence_score"] = np.asarray(Xn[:, cols].mean(axis=1)).ravel()

    fig, ax = plt.subplots(figsize=(5, 5))
    data = [adata.obs.loc[adata.obs.age_group == g, "senescence_score"].values
            for g in ["young", "old"]]
    ax.violinplot(data, showmedians=True)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["young", "old"])
    ax.set_ylabel("senescence score ({} markers, log-norm)".format(len(present)))
    ax.set_title("Senescence marker score by age")
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "senescence_violin.png"), dpi=150)
    plt.close(fig)

    by = adata.obs.groupby(["_celltype", "age_group"])["senescence_score"].mean().reset_index()
    by.to_csv(os.path.join(data_dir, "senescence_by_celltype.csv"), index=False)
    print("[senescence] used markers: {}".format(", ".join(present)))
    return present


def umap(adata, results_dir):
    if "X_umap" not in adata.obsm:
        print("[umap] no precomputed X_umap, skipping (compute with scanpy.tl.umap if needed)")
        return False
    emb = np.asarray(adata.obsm["X_umap"])
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    cts = adata.obs["_celltype"].astype("category")
    for name in cts.cat.categories:
        m = (cts == name).to_numpy()
        axs[0].scatter(emb[m, 0], emb[m, 1], s=2, alpha=0.5, label=str(name))
    axs[0].set_title("UMAP by cell type")
    axs[0].axis("off")
    if len(cts.cat.categories) <= 14:
        axs[0].legend(markerscale=4, fontsize=6, loc="best")
    for grp in ["young", "old"]:
        m = (adata.obs.age_group == grp).to_numpy()
        axs[1].scatter(emb[m, 0], emb[m, 1], s=2, alpha=0.4, color=COLORS[grp], label=grp)
    axs[1].set_title("UMAP by age")
    axs[1].axis("off")
    axs[1].legend(markerscale=4)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "umap.png"), dpi=150)
    plt.close(fig)
    return True


def volcano(res, cell_type, results_dir):
    r = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    r["nl"] = -np.log10(r["padj"].clip(lower=1e-300))
    sig = r["padj"] < 0.05
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(r.loc[~sig, "log2FoldChange"], r.loc[~sig, "nl"], s=8, alpha=0.4,
               color="#bdbdbd", edgecolor="none")
    ax.scatter(r.loc[sig, "log2FoldChange"], r.loc[sig, "nl"], s=10, alpha=0.7,
               color="#cb181d", edgecolor="none")
    ax.axhline(-np.log10(0.05), color="black", lw=0.8, ls="--")
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("log2 fold change (old vs young)")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title("Pseudobulk DE: {}".format(cell_type))
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "pseudobulk_volcano.png"), dpi=150)
    plt.close(fig)
