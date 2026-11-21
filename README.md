# Breast Cancer Gene Expression (TCGA-BRCA)

Can gene-expression profiles distinguish breast **tumor** tissue from **normal** breast tissue, and which genes give it away?

This project analyzes RNA-seq data from **TCGA-BRCA** (The Cancer Genome Atlas, Breast Invasive Carcinoma), downloaded from the NCI [Genomic Data Commons](https://portal.gdc.cancer.gov/). It is a learning and portfolio project built step by step: data acquisition → quality control → differential expression → machine learning → biological interpretation.

> 🚧 **Status:** work in progress. Project setup is done; data acquisition is next.

## Data

- **Source:** GDC, project `TCGA-BRCA`, RNA-seq, *STAR - Counts* files (one TSV per sample, gene model GENCODE v36)
- **Samples:** ~1,110 Primary Tumor and ~113 Solid Tissue Normal
- **Metadata:** the GDC sample sheet (maps each file to its patient and sample type)

The data is not stored in this repository. See [Reproducing the analysis](#reproducing-the-analysis).

## Project structure

```
├── data/
│   ├── raw/          # GDC downloads (not in Git)
│   └── processed/    # merged matrices (not in Git)
├── notebooks/        # exploration and figures
├── src/              # reusable pipeline code
├── results/
│   └── figures/
└── requirements.txt
```

## Reproducing the analysis

Requires Python 3.14.

```bash
python -m venv .venv
.venv\Scripts\activate        # Windows (macOS/Linux: source .venv/bin/activate)
pip install -r requirements.txt
```

Data download and pipeline steps will be added as each stage is built.

## Methods (planned)

1. Build gene-by-sample matrices of raw counts and TPM
2. Quality control and PCA on log2(TPM + 1)
3. Differential expression on raw counts, with Benjamini-Hochberg FDR correction
4. Tumor vs. normal classification, with splits by patient to avoid leakage, and AUROC / PR-AUC / balanced accuracy reported because of class imbalance
5. Biological interpretation of the top genes

## Results

*Coming soon.*
