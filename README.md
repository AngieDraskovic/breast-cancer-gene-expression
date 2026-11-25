# Breast Cancer Gene Expression (TCGA-BRCA)

Can gene-expression profiles distinguish breast **tumor** tissue from **normal** breast tissue, and which genes give it away?

This project analyzes RNA-seq data from **TCGA-BRCA** (The Cancer Genome Atlas, Breast Invasive Carcinoma), downloaded from the NCI [Genomic Data Commons](https://portal.gdc.cancer.gov/). It is a learning and portfolio project built step by step: data acquisition → quality control → differential expression → machine learning → biological interpretation.

> 🚧 **Status:** work in progress. Data acquisition and matrix building are done; quality control is next.

## Data

- **Source:** GDC, project `TCGA-BRCA`, RNA-seq, *STAR - Counts* files (one TSV per sample, gene model GENCODE v36)
- **Files:** 1,231 open-access files from 1,095 patients: 1,111 Primary Tumor, 113 Solid Tissue Normal, 7 Metastatic
- **Metadata:** [`metadata/sample_sheet.tsv`](metadata/sample_sheet.tsv), built from the GDC API, maps each file to its patient and sample type

The data is not stored in this repository. See [Reproducing the analysis](#reproducing-the-analysis).

## Project structure

```
├── data/
│   ├── raw/          # GDC downloads (not in Git)
│   └── processed/    # merged matrices (not in Git)
├── metadata/         # sample sheet: which GDC files were analyzed
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

Then, from the project root:

```bash
python src/download_gdc.py      # sample sheet + 1,231 count files (~5 GB) into data/raw/
python src/build_matrices.py    # counts/TPM matrices + sample metadata into data/processed/
```

## Methods

### Data preparation

- For each sample, raw read counts (`unstranded`) and TPM (`tpm_unstranded`) were extracted into two separate gene-by-sample matrices (60,660 genes × 1,226 samples), after removing STAR's four summary rows (`N_unmapped`, etc.).
- 5 tumor samples were sequenced twice; for each, we kept the run with more reads assigned to genes, as it gives the more precise measurement. Every dropped file is logged.
- Both matrices are built in a single pass over the files: each file is read once, and the duplicate-sample rule is applied in one place, so the counts and TPM matrices are guaranteed to come from the same files.
- The 7 metastatic samples are kept in the matrices but excluded from the tumor-vs-normal analyses.
- Automated sanity checks confirm the matrix shapes, matching sample order across matrices and metadata, expected group sizes, no missing or negative values, and TPM columns summing to ~1,000,000.

### Planned

1. Quality control and PCA
3. Differential expression on raw counts, with Benjamini-Hochberg FDR correction
4. Tumor vs. normal classification, with splits by patient to avoid leakage, and AUROC / PR-AUC / balanced accuracy reported because of class imbalance
5. Biological interpretation of the top genes

## Results

*Coming soon.*
