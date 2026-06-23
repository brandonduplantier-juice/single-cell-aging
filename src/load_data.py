"""
Load a Tabula Muris Senis tissue h5ad and label cells young vs old.

The sandbox cannot download the data. Get one tissue file from Tabula Muris Senis
(CZ CELLxGENE at cellxgene.cziscience.com, or the TMS figshare collection) and place
the .h5ad in data/. The pipeline reads whichever single .h5ad is in data/.

Column names differ between TMS releases, so this auto-detects the age, cell-type,
and mouse columns and prints which it used. If detection fails, pass them explicitly.
"""

import glob
import os
import re

import numpy as np
import pandas as pd
import anndata as ad

AGE_COLS = ["age", "development_stage", "Age", "mouse.age"]
CELLTYPE_COLS = ["cell_type", "cell_ontology_class", "free_annotation",
                 "celltype", "cell_ontology_class_reannotated"]
MOUSE_COLS = ["mouse.id", "donor_id", "mouse_id", "donor", "mouse"]
YOUNG_MAX_MONTHS = 6
OLD_MIN_MONTHS = 18


def _pick(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def _months(val):
    s = str(val).lower()
    m = re.search(r"(\d+)\s*(?:m|month)", s)
    if m:
        return int(m.group(1))
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else None


def find_h5ad(data_dir):
    files = sorted(glob.glob(os.path.join(data_dir, "*.h5ad")))
    if not files:
        raise FileNotFoundError(
            "no .h5ad in {}. Download a Tabula Muris Senis tissue file there "
            "(cellxgene.cziscience.com or the TMS figshare).".format(data_dir))
    return files[0]


def load(path, age_col=None, celltype_col=None, mouse_col=None):
    adata = ad.read_h5ad(path)
    cols = list(adata.obs.columns)
    age_col = age_col or _pick(cols, AGE_COLS)
    celltype_col = celltype_col or _pick(cols, CELLTYPE_COLS)
    mouse_col = mouse_col or _pick(cols, MOUSE_COLS)
    missing = [n for n, v in [("age", age_col), ("cell_type", celltype_col),
                              ("mouse", mouse_col)] if v is None]
    if missing:
        raise KeyError("could not auto-detect columns {}. Available: {}. "
                       "Pass them explicitly to load().".format(missing, cols))

    months = pd.to_numeric(adata.obs[age_col].astype(str).map(_months), errors="coerce")
    grp = np.where(months <= YOUNG_MAX_MONTHS, "young",
                   np.where(months >= OLD_MIN_MONTHS, "old", "mid"))
    adata.obs["age_group"] = grp
    adata.obs["_celltype"] = adata.obs[celltype_col].astype(str)
    adata.obs["_mouse"] = adata.obs[mouse_col].astype(str)
    adata = adata[adata.obs["age_group"].isin(["young", "old"])].copy()
    print("[load] {} | age={} celltype={} mouse={} | young {} / old {} cells".format(
        os.path.basename(path), age_col, celltype_col, mouse_col,
        int((adata.obs.age_group == "young").sum()),
        int((adata.obs.age_group == "old").sum())))
    return adata, {"age": age_col, "celltype": celltype_col, "mouse": mouse_col}
