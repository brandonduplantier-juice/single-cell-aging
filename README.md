# Single-Cell Aging Pathway Analysis (Tabula Muris Senis)

Finds which cell types shift and which genes change with age in mouse tissue, from
single-cell data, analyzed the statistically correct way: pseudobulk, with mice as
the replicates.

Version 1.0.0

## What it shows (after you run it on a tissue)

For one mouse tissue, young versus old, it produces three things: how the mix of cell
types shifts with age, a senescence (cellular-aging) score per cell type, and the
genes that get louder or quieter with age within each cell type. The numbers come
from your run; this build cannot download the atlas from a sandbox.

Figures: results/umap.png, results/composition.png, results/senescence_violin.png,
results/pseudobulk_volcano.png

## In plain English

Here is the idea, and the one technical choice that makes or breaks it, with terms
defined as they appear.

**What single-cell data is.** Older methods measured the average gene activity across
a whole lump of tissue. *Single-cell* sequencing instead measures gene activity in
each individual cell, thousands of them, so you can see that a tissue is a mix of
different *cell types* (immune cells, structural cells, and so on), each behaving
differently.

**What this project asks.** In a given mouse tissue, what changes between young and
old animals? Two things: does the *mix* of cell types shift (more of some, fewer of
others), and within a cell type, which genes change their activity.

**The trap that ruins most amateur versions.** Because there are thousands of cells
per mouse, it is tempting to treat each cell as a separate data point. That is
*pseudoreplication*: cells from the same animal are not independent, so counting them
as independent produces wildly overconfident, basically fake, statistics. A reviewer
checks for this first.

**The correct fix, and what this repo does.** *Pseudobulk*: add up each animal's cells
within a cell type into one profile per mouse, then compare across mice. Now the unit
of comparison is the animal, which is the honest replicate. Doing pseudobulk, and
being able to explain why, is the entire difference between a credible single-cell
aging result and a misleading one.

**The extra aging readout.** It also scores a panel of *senescence* markers.
Senescence is a state where cells stop dividing but linger and cause low-grade
inflammation, a hallmark of aging. A score that rises with age is the expected,
on-theme signal.

What this is **not**: it is mouse, not human; it analyzes pre-labeled public data
rather than generating it; and cell-type proportions are relative (if one type rises,
others mechanically fall), so a proportion shift is not by itself proof a cell type's
absolute number changed.

## How it works (technical)

Loads one tissue .h5ad, plots the precomputed UMAP by cell type and age, computes
cell-type composition per mouse, scores senescence and SASP markers per cell with
scanpy, and runs pydeseq2 pseudobulk differential expression (old vs young) within
each cell type that has at least two mice per group. It uses the atlas's existing cell
labels rather than re-clustering, to spend effort on the analysis, not on redoing the
annotation.

## What is in this repo

    run_all.py            the full pipeline (needs one tissue .h5ad in data/)
    src/                  loader, pseudobulk aggregation, analysis
    app.py                light dashboard that reads small precomputed summaries
    requirements-app.txt  minimal deps so the dashboard deploys without the big file

## Data

Tabula Muris Senis, a mouse aging cell atlas. Download one tissue .h5ad from CZ
CELLxGENE (cellxgene.cziscience.com) or the TMS figshare collection and place it in
data/. Young defaults to 3 months or less, old to 18 months or more. The loader
auto-detects the age, cell-type, and mouse columns across TMS releases and prints
which it used.

## Run the pipeline (Windows)

    python -m venv .venv
    .venv\Scripts\python.exe -m pip install -r requirements.txt
    (put one TMS tissue .h5ad in data\)
    .venv\Scripts\python.exe run_all.py

## Explore it (after a run)

    pip install -r requirements-app.txt
    python -m streamlit run app.py

## Limitations

Mouse, not human. Pre-annotated public data, not generated here. Composition shifts
are relative, not absolute counts. Cross-condition single-cell DE is valid only at the
pseudobulk level, which is why it is built that way. Findings are framed as consistent
with known aging biology, not as discoveries.

## Glossary

- **Single-cell sequencing**: measuring gene activity in each individual cell, not the
  tissue average.
- **Cell type**: a category of cell (for example immune or structural) with a distinct
  activity profile.
- **Pseudoreplication**: wrongly treating non-independent units (cells from one animal)
  as independent, which fakes significance.
- **Pseudobulk**: summing each animal's cells into one profile per animal, the correct
  unit for comparison.
- **Composition shift**: a change in the proportions of cell types between groups.
- **Senescence / SASP**: a non-dividing, inflammation-causing cell state tied to aging.
- **UMAP**: a 2D map that places similar cells near each other for visualization.
- **h5ad**: the standard file format that stores a single-cell dataset.

## Citation and disclaimer

Data: Tabula Muris Senis (The Tabula Muris Consortium, Nature 2020); raw GSE132042.
Method: pseudobulk DE via pydeseq2. Research and portfolio project.
