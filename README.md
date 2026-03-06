# Batch Fine-mapping and Colocalization Pipeline



This repository contains an automated R pipeline for processing Genome-Wide Association Study (GWAS) summary statistics. It performs greedy LD clumping to identify independent loci, runs fine-mapping using **FINEMAP**, and performs colocalization analysis using the **coloc** R package.

## 🚀 Features
* **Automated Locus Partitioning:** Identifies genome-wide significant lead SNPs and extracts a configurable genomic window (default +/- 1.5MB).
* **FINEMAP Integration:** Automatically generates `.z` and `.master` files and executes the FINEMAP binary for each locus.
* **Colocalization:** Integrates with QTL data to compute the posterior probability (H4) of a shared causal variant.
* **Aggregated Outputs:** Combines results into clean `.csv`, `.xlsx`, and `.txt` files for easy downstream analysis.

## 📋 Prerequisites

### 1. Software Requirements
* **R (>= 3.6.0)**
* **FINEMAP binary:** Ensure FINEMAP is installed and executable on your system. 

### 2. R Packages
The script will attempt to install missing packages automatically, but requires:
`argparse`, `data.table`, `dplyr`, `ggplot2`, `plyr`, `devtools`, `stringr`, `tibble`, `coloc`, `writexl`

### 3. LD Matrices
You must have pre-computed LD correlation matrices for your loci. By default, the script looks in `./ld_matrices/` for files named `<rsid>_ld.txt` (e.g., `rs12345_ld.txt`).

## 💻 Usage

Run the pipeline from the command line using `Rscript`:

```bash
Rscript run_pipeline.R \
  --z_file "data/my_gwas_data.csv" \
  --qtl_file "data/my_eqtl_data.csv" \
  --study_id "HeartDisease_GWAS" \
  --pip_threshold 0.5 \
  --finemap_exe "/path/to/finemap"
