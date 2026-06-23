"""
Pseudobulk differential expression, the correct way to test single-cell DE.

Aggregates raw counts to one profile per (cell type, mouse), so the statistical
replicates are animals, not cells. Per-cell DE treats thousands of cells from one
mouse as independent samples, which is pseudoreplication and produces tiny, fake
p-values. Pseudobulk avoids that. Runs pydeseq2 (old vs young) within each cell type
that has enough mice per group.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_MICE_PER_GROUP = 2
MIN_CELLS_PER_PSEUDOSAMPLE = 10
MIN_GENE_TOTAL = 10


def get_counts(adata):
    """Return a raw-count matrix and the matching gene names."""
    if "counts" in adata.layers:
        return adata.layers["counts"], list(adata.var_names)
    if adata.raw is not None:
        return adata.raw.X, list(adata.raw.var_names)
    return adata.X, list(adata.var_names)


def make_pseudobulk(adata):
    X, var_names = get_counts(adata)
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    idx = pd.DataFrame({
        "cell_type": adata.obs["_celltype"].to_numpy(),
        "mouse": adata.obs["_mouse"].to_numpy(),
        "age_group": adata.obs["age_group"].to_numpy(),
    })
    rows, meta = [], []
    for (ct, mouse, age), sub in idx.groupby(["cell_type", "mouse", "age_group"]):
        ii = sub.index.to_numpy()
        if len(ii) < MIN_CELLS_PER_PSEUDOSAMPLE:
            continue
        rows.append(np.asarray(X[ii].sum(axis=0)).ravel())
        meta.append((ct, mouse, age, len(ii)))
    pb = pd.DataFrame(rows, columns=var_names)
    md = pd.DataFrame(meta, columns=["cell_type", "mouse", "age_group", "n_cells"])
    return pb, md


def run_de_per_celltype(pb, md):
    """Old-vs-young pydeseq2 within each cell type with enough mice. Returns dict."""
    results = {}
    for ct in sorted(md["cell_type"].unique()):
        m = (md["cell_type"] == ct).to_numpy()
        sub_md = md[m].reset_index(drop=True)
        n_young = int((sub_md.age_group == "young").sum())
        n_old = int((sub_md.age_group == "old").sum())
        if n_young < MIN_MICE_PER_GROUP or n_old < MIN_MICE_PER_GROUP:
            continue
        counts = pb[m].reset_index(drop=True).round().astype(int)
        counts = counts.loc[:, counts.sum(axis=0) >= MIN_GENE_TOTAL]
        if counts.shape[1] < 10:
            continue
        design = sub_md[["age_group"]].copy()
        dds = DeseqDataSet(counts=counts, metadata=design,
                           design_factors="age_group", quiet=True)
        dds.deseq2()
        stat = DeseqStats(dds, contrast=["age_group", "old", "young"], quiet=True)
        stat.summary()
        res = stat.results_df.copy()
        res["symbol"] = res.index
        results[ct] = res.sort_values("padj")
        print("  [de] {}: {} young / {} old mice, {} genes tested".format(
            ct, n_young, n_old, counts.shape[1]))
    return results
