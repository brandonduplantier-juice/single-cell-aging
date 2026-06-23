# Single-Cell Aging Pathway Analysis (Tabula Muris Senis)

Which cell types change, and which genes shift with age, in mouse tissue single-cell
data, analyzed the statistically correct way: pseudobulk, with mice as the replicates.

Version 1.0.0

## Why pseudobulk (the part that matters)

Single-cell data has thousands of cells per animal. Testing differential expression
per cell treats those thousands as independent samples, which is pseudoreplication and
produces falsely tiny p-values. This pipeline aggregates counts to one profile per
mouse per cell type and tests across mice, so the replicate is the animal. That is the
difference between a credible single-cell aging result and a misleading one, and it is
the first thing a reviewer checks.

## Result

Run it on a tissue and results/metrics.json plus the figures summarize three things:
how cell-type composition shifts young vs old, a senescence-marker score that rises
with age, and pseudobulk differential expression per cell type. The numbers come from
your run; this sandbox cannot download the atlas.

Figures: results/umap.png, results/composition.png, results/senescence_violin.png,
results/pseudobulk_volcano.png

## Explore it

The dashboard reads small precomputed summary tables, so it deploys without the large
h5ad:

    pip install -r requirements-app.txt
    python -m streamlit run app.py

It shows cell-type composition shifts, senescence by cell type, and the genes that get
louder or quieter with age in each cell type, all in plain language.

## Data

Tabula Muris Senis, a mouse aging cell atlas. Download one tissue .h5ad from CZ
CELLxGENE (cellxgene.cziscience.com) or the TMS figshare collection and place it in
data/. By default young is 3 months or less and old is 18 months or more. The loader
auto-detects the age, cell-type, and mouse columns across TMS releases and prints which
it used.

## Run the pipeline (Windows)

    python -m venv .venv
    .venv\Scripts\python.exe -m pip install -r requirements.txt
    (put one TMS tissue .h5ad in data\)
    .venv\Scripts\python.exe run_all.py

## Method

Loads the tissue, plots the precomputed UMAP by cell type and age, computes cell-type
composition per mouse (so differences reflect animals, not single noisy cells), scores
a panel of senescence and SASP markers per cell, and runs pydeseq2 pseudobulk DE old vs
young within each cell type that has at least two mice per group.

## Outputs

results/: umap.png, composition.png, senescence_violin.png, pseudobulk_volcano.png,
metrics.json. data/: composition.csv, senescence_by_celltype.csv, pseudobulk_top.csv,
and de/<cell_type>.csv per cell type.

## Three things to confirm before trusting the numbers

1. Raw counts. Pseudobulk needs integer counts, not normalized values. The loader
   prefers layers['counts'] or .raw; confirm your file has raw counts.
2. Column detection. Check the age, cell-type, and mouse columns printed at load. If
   wrong, pass them to load_data.load() explicitly.
3. Age thresholds. The default young (<=6 months) and old (>=18 months) cutoffs should
   match the ages actually sampled in your tissue.

## Limitations

Cell-type labels come from the source atlas. Cell types with too few mice per group are
skipped. The senescence score is a simple marker-panel mean, not a validated assay. This
is mouse data, so human relevance is indirect. UMAP uses the atlas embedding. Research
and education only.

## Citation

Data: The Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes
ageing tissues in the mouse. Nature, 2020.
Method: Love MI, Huber W, Anders S. DESeq2. Genome Biology, 2014, applied to pseudobulk
profiles.
