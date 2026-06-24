"""
Single-cell aging explorer.

A plain-language dashboard over precomputed summaries of a Tabula Muris Senis tissue:
how cell-type makeup shifts with age, how a senescence score changes, and which genes
shift with age in each cell type. Reads only small CSVs, so it runs without the large
h5ad. Mouse data, research and education only, not medical advice.

Run locally:
  pip install -r requirements-app.txt
  python -m streamlit run app.py
"""

import json
import os

import pandas as pd
import plotly.express as px
import streamlit as st

APP_VERSION = "1.0.0"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RESULTS = os.path.join(HERE, "results")
COLORS = {"young": "#2c7fb8", "old": "#d95f0e"}


@st.cache_data
def load():
    comp = pd.read_csv(os.path.join(DATA, "composition.csv"))
    sen = _maybe(os.path.join(DATA, "senescence_by_celltype.csv"))
    top = _maybe(os.path.join(DATA, "pseudobulk_top.csv"))
    meta = {}
    mp = os.path.join(RESULTS, "metrics.json")
    if os.path.exists(mp):
        with open(mp) as fh:
            meta = json.load(fh)
    return comp, sen, top, meta


def _maybe(path):
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


if not os.path.exists(os.path.join(DATA, "composition.csv")):
    st.error("Summaries not found. Run: python run_all.py")
    st.stop()

comp, sen, top, meta = load()

st.title("How tissue ages, cell by cell")
st.caption("Cell-type and gene-activity changes between young and old mice, from the "
           "Tabula Muris Senis atlas. Research and education only.")

with st.expander("New here? How to read this", expanded=False):
    st.markdown(
        "- Single-cell data measures gene activity in thousands of individual cells, "
        "each labeled with its **cell type**.\n"
        "- **Composition** is the mix of cell types. Aging can shift the mix (for "
        "example more immune cells, fewer stem cells).\n"
        "- **Senescence** is a worn-out cell state that builds up with age. The score "
        "here averages a panel of known senescence markers.\n"
        "- **Aging genes** are genes that get louder or quieter with age within a cell "
        "type, found by comparing old vs young mice.\n"
        "- These comparisons are done per **mouse**, not per cell, which is the honest "
        "way to avoid counting the same animal thousands of times."
    )

with st.expander("Key terms"):
    st.markdown(
        "- **Single-cell RNA-seq**: measuring gene activity in thousands of individual "
        "cells separately, instead of averaging a whole tissue.\n"
        "- **Cell type**: the category a cell belongs to, like immune or stem cell.\n"
        "- **Composition**: the mix of cell types. Aging can shift the mix.\n"
        "- **Aging genes**: genes that get louder or quieter with age within one cell "
        "type.\n"
        "- **Senescence**: a worn-out cell state that builds up with age. The score "
        "averages a panel of known markers.\n"
        "- **Pseudobulk**: grouping each mouse's cells together before testing, so the "
        "statistics count animals, not cells. The honest way to avoid counting one "
        "animal thousands of times.\n"
        "- **FDR**: the false-discovery-rate cutoff (0.05 here) that keeps false alarms "
        "low when testing many genes at once."
    )

if meta:
    c1, c2, c3 = st.columns(3)
    c1.metric("Cells analyzed", "{:,}".format(meta.get("n_cells", 0)))
    c2.metric("Cell types", meta.get("n_cell_types", comp["_celltype"].nunique()))
    total_sig = sum(meta.get("n_sig_by_celltype", {}).values())
    c3.metric("Age-changed genes", total_sig,
              help="Total genes significantly shifted with age across cell types (pseudobulk, FDR 0.05).")

st.subheader("Which cell types grow or shrink with age")
st.caption("Each bar is a cell type's share of all cells. Compare young (blue) and old "
           "(orange): a longer orange bar means that cell type becomes more common with age.")
fig = px.bar(comp, x="frac", y="_celltype", color="age_group", barmode="group",
             color_discrete_map=COLORS, orientation="h",
             labels={"frac": "fraction of cells", "_celltype": "", "age_group": "age"})
fig.update_layout(height=max(320, 26 * comp["_celltype"].nunique()))
st.plotly_chart(fig, use_container_width=True)

if len(sen):
    st.subheader("Where senescence builds up")
    st.caption("Higher means more of the worn-out cell signature. Compare young vs old "
               "within each cell type.")
    sfig = px.bar(sen, x="senescence_score", y="_celltype", color="age_group",
                  barmode="group", color_discrete_map=COLORS, orientation="h",
                  labels={"senescence_score": "senescence score", "_celltype": "", "age_group": "age"})
    sfig.update_layout(height=max(320, 26 * sen["_celltype"].nunique()))
    st.plotly_chart(sfig, use_container_width=True)

st.subheader("Aging genes in one cell type")
cell_types = sorted(comp["_celltype"].unique())
pick = st.selectbox("Cell type", cell_types)

cc = comp[comp["_celltype"] == pick]
young_f = float(cc.loc[cc.age_group == "young", "frac"].sum())
old_f = float(cc.loc[cc.age_group == "old", "frac"].sum())
m1, m2 = st.columns(2)
m1.metric("Share of cells, young", "{:.1%}".format(young_f))
m2.metric("Share of cells, old", "{:.1%}".format(old_f),
          delta="{:+.1%}".format(old_f - young_f))

if len(top):
    t = top[top["cell_type"] == pick].copy()
    if len(t):
        t["effect"] = t["log2FoldChange"].apply(
            lambda x: "about {:.1f}x {}".format(2 ** abs(x), "louder" if x > 0 else "quieter"))
        up = t[t["log2FoldChange"] > 0].sort_values("padj")
        down = t[t["log2FoldChange"] < 0].sort_values("padj")
        a, b = st.columns(2)
        a.markdown("**Louder in old age**")
        a.dataframe(up[["symbol", "effect", "padj"]], hide_index=True, use_container_width=True)
        b.markdown("**Quieter in old age**")
        b.dataframe(down[["symbol", "effect", "padj"]], hide_index=True, use_container_width=True)
    else:
        st.info("No significant aging genes for this cell type (or too few mice to test it).")

st.caption("Summaries precomputed from a mouse single-cell atlas. Differential "
           "expression is pseudobulk (per mouse), not per cell. Research and education "
           "only, not medical advice. App v{}.".format(APP_VERSION))
