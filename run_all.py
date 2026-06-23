"""
Run the single-cell aging pipeline: load, composition, senescence, UMAP, pseudobulk DE.

Put one Tabula Muris Senis tissue .h5ad in data/ first (see src/load_data.py), then:
  python run_all.py
"""

import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))

import load_data      # noqa: E402
import pseudobulk     # noqa: E402
import analysis       # noqa: E402

__version__ = "1.0.0"
DATA = os.path.join(HERE, "data")
RESULTS = os.path.join(HERE, "results")


def main():
    print("[run_all] single-cell-aging v{}".format(__version__))
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(os.path.join(DATA, "de"), exist_ok=True)

    path = load_data.find_h5ad(DATA)
    adata, colmap = load_data.load(path)

    analysis.umap(adata, RESULTS)
    analysis.composition(adata, RESULTS, DATA)
    analysis.senescence(adata, RESULTS, DATA)

    pb, md = pseudobulk.make_pseudobulk(adata)
    res = pseudobulk.run_de_per_celltype(pb, md)

    top_rows = []
    for ct, r in res.items():
        safe = ct.replace("/", "_").replace(" ", "_")
        r.to_csv(os.path.join(DATA, "de", safe + ".csv"))
        sig = r[r["padj"] < 0.05]
        for _, row in sig.head(20).iterrows():
            top_rows.append({"cell_type": ct, "symbol": row["symbol"],
                             "log2FoldChange": row["log2FoldChange"], "padj": row["padj"]})
    pd.DataFrame(top_rows).to_csv(os.path.join(DATA, "pseudobulk_top.csv"), index=False)

    if res:
        focal = max(res, key=lambda k: int((res[k]["padj"] < 0.05).sum()))
        analysis.volcano(res[focal], focal, RESULTS)

    metrics = {
        "n_cells": int(adata.n_obs),
        "n_cell_types": int(adata.obs["_celltype"].nunique()),
        "de_cell_types": list(res.keys()),
        "n_sig_by_celltype": {ct: int((r["padj"] < 0.05).sum()) for ct, r in res.items()},
    }
    with open(os.path.join(RESULTS, "metrics.json"), "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, indent=2)
    print("[run_all] done. Summaries in data/, figures in results/. "
          "Launch: python -m streamlit run app.py")


if __name__ == "__main__":
    main()
