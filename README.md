# FaSTLMM_Runner

Scripts to run Genome-Wide Association Studies (GWAS) using FaST-LMM and perform permutation tests to detect QTL.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Workflow](#workflow)
- [Scripts Reference](#scripts-reference)
  - [simpleGWAS\_VV\_Template.sh](#simplegwas_vv_templatesh)
  - [runParallel.sh](#runparallelsh)
  - [src/runFaSTLMM.py](#srcrunfastlmmpy)
  - [src/rank\_based\_Inverse\_Normal\_Transformation.R](#srcrank_based_inverse_normal_transformationr)
  - [src/calc\_GIF.R](#srccalc_gifr)
  - [src/qqman\_script.R](#srcqqman_scriptr)
  - [src/addLinkageGroups.py](#srcaddlinkagegroupspy)
- [Input File Formats](#input-file-formats)
- [Output Files](#output-files)
- [HPC Configuration](#hpc-configuration)
- [References](#references)

---

## Overview

FaSTLMM_Runner orchestrates a complete GWAS pipeline around the [FaST-LMM](https://github.com/fastlmm/FaST-LMM) software (tested with v0.6.4). Starting from raw phenotype and genotype data, the pipeline:

1. Normalises phenotypes using a rank-based inverse normal transformation (INT).
2. Runs FaST-LMM association testing on each normalised phenotype.
3. Estimates an empirical significance threshold via permutation testing (5th percentile of best per-permutation p-values).
4. Computes the Genomic Inflation Factor (λ) to assess p-value calibration.
5. Generates Manhattan and QQ plots.
6. Annotates significant SNPs with linkage group information using a graph-based LD approach.

The pipeline is designed for an HPC cluster using the CEA/TGCC `ccc_msub`/`ccc_mprun` scheduler, but the core Python and R scripts are scheduler-agnostic and can be adapted to SLURM, PBS, or local execution.

---

## Repository Structure

```
FaSTLMM_Runner/
├── simpleGWAS_VV_Template.sh                       # Main HPC job script (single run)
├── runParallel.sh                                  # Launcher: splits phenotypes across N parallel jobs
└── src/
    ├── runFaSTLMM.py                               # FaST-LMM association + permutation test (v4.0)
    ├── rank_based_Inverse_Normal_Transformation.R  # Phenotype normalisation
    ├── calc_GIF.R                                  # Genomic Inflation Factor calculation
    ├── qqman_script.R                              # Manhattan and QQ plot generation
    └── addLinkageGroups.py                         # LD-based linkage group annotation
```

---

## Requirements

### System
- Bash
- [PLINK 1.9](https://www.cog-genomics.org/plink/) — required by `addLinkageGroups.py` for pairwise LD computation
- HPC environment with `ccc_msub`/`ccc_mprun` (CEA TGCC), or adapt submission commands for your scheduler

### Python (≥ 3.8)

| Package | Purpose |
|---------|---------|
| `fastlmm` ≥ 0.6.4 | GWAS association testing |
| `pysnptools` | Genotype/phenotype data I/O |
| `numpy` | Numerical operations |
| `pandas` | Data manipulation |
| `networkx` | Graph-based linkage group clustering |

```bash
pip install fastlmm pysnptools numpy pandas networkx
```

### R

| Package | Purpose |
|---------|---------|
| `qqman` | Manhattan and QQ plots |
| `stringr` | String manipulation (normalisation script) |
| `calibrate` | Optional, loaded in `qqman_script.R` |

---

## Quick Start

1. **Clone the repository:**
   ```bash
   git clone https://github.com/VLoegler/FaSTLMM_Runner.git
   cd FaSTLMM_Runner
   ```

2. **Prepare input data** (see [Input File Formats](#input-file-formats)):
   - Genotypes in PLINK binary format (`.bed`, `.bim`, `.fam`)
   - One `.phen` file per trait in a shared directory
   - Optional: kinship matrix (PLINK format) and/or covariance matrix

3. **Edit `simpleGWAS_VV_Template.sh`** — fill in the variables block:
   ```bash
   WORKDIR=/path/to/workdir
   GENO=/path/to/genotypes/genotypes.plink
   PHENODIR=/path/to/phenotypes
   KINSHIP=""      # or /path/to/kinship.plink
   COVAR=""        # or /path/to/covariance.plink
   NB_PERM=100
   SCRIPT_DIR=/path/to/src
   ```

4. **Submit a single job:**
   ```bash
   ccc_msub simpleGWAS_VV_Template.sh
   ```

5. **Or parallelise** across multiple jobs (edit `runParallel.sh` first):
   ```bash
   bash runParallel.sh
   ```

---

## Workflow

```
PHENODIR/*.phen
      │
      ▼  [rank_based_Inverse_Normal_Transformation.R]
      │     Blom INT (c=3/8); NAs removed from output
      ▼
PHENODIR/*.norm.phen
      │
      ▼  [runFaSTLMM.py]  ◄── GENO, optional KINSHIP, COVAR
      │     Real association via fastlmm.single_snp
      │     NB_PERM permutations, batched ≤ 101 per call
      │     Threshold = 5th percentile of best permutation p-values
      │
      ├─► <COND>.first_assoc.txt.gz    (all SNP results)
      ├─► <COND>.threshold.txt         (empirical p-value threshold)
      └─► <COND>.signif_snps.txt       (SNPs below threshold)
              │
              ├── [calc_GIF.R]          → <COND>.first_assoc.txt.gz.lgc  (λ value)
              ├── [qqman_script.R]      → Manhattan (.eps/.jpg) + QQ (.jpg) plots
              └── [addLinkageGroups.py] → <COND>.signif_snps.LinkageGroups.txt
```

---

## Scripts Reference

### `simpleGWAS_VV_Template.sh`

The main pipeline script, submitted as a single HPC batch job. Runs all three steps (normalisation → GWAS → post-processing) sequentially across all phenotypes in `PHENODIR`.

**HPC directives (MSUB):**

| Directive | Value | Meaning |
|-----------|-------|---------|
| `-r simple_GWAS` | — | Job name |
| `-n 1` | 1 node | |
| `-c 8` | 8 cores | Forwarded to `runFaSTLMM.py -t` |
| `-Q long` | long queue | |
| `-T 259200` | 72 h | Walltime in seconds |
| `-q milan` | milan | Target partition |

**Variables to configure:**

| Variable | Description |
|----------|-------------|
| `WORKDIR` | Root working directory |
| `GENO` | PLINK genotype prefix (no extension) |
| `PHENODIR` | Directory containing `.phen` / `.norm.phen` files |
| `KINSHIP` | PLINK kinship prefix, or `""` to use genotypes as GRM |
| `COVAR` | Covariance matrix file path, or `""` for none |
| `NB_PERM` | Number of permutations for threshold estimation |
| `SCRIPT_DIR` | Path to the `src/` directory |

**Step 1 — Normalisation:** For each `.phen` file in `PHENODIR` that does not already have a `.norm.phen` counterpart, runs `rank_based_Inverse_Normal_Transformation.R`.

**Step 2 — GWAS:** For each `.norm.phen`, calls `runFaSTLMM.py` via `ccc_mprun`, with `-k` and/or `-c` flags set conditionally based on `KINSHIP` and `COVAR`.

**Step 3 — Post-processing:** For each phenotype, reads the condition name from the first line of the `.norm.phen` file and runs `qqman_script.R`, `calc_GIF.R`, and `addLinkageGroups.py`.

---

### `runParallel.sh`

Distributes phenotypes across `NB_PARA` parallel HPC jobs to reduce overall wall-clock time.

**How it works:**
1. Creates `NB_PARA` subdirectories (`phenotypes_Para1` … `phenotypes_ParaN`) under `WORKDIR`.
2. Round-robin copies `.phen` files (and companions) from `PHENODIR` into those directories.
3. For each parallel slot, generates a customised copy of `simpleGWAS_VV_Template.sh` by injecting variable values with `sed`, then submits it with `ccc_msub`. Results go to separate `*_fastlmm_results_ParaN` directories.

**Variables to configure:**

| Variable | Description |
|----------|-------------|
| `NB_PARA` | Number of parallel jobs (default: `5`) |
| `WORKDIR` | Root working directory |
| `GENO` | PLINK genotype prefix |
| `PHENODIR` | Source phenotype directory |
| `KINSHIP` | PLINK kinship prefix, or `""` |
| `COVAR` | Covariance matrix file, or `""` |
| `NB_PERM` | Number of permutations |

> **Note:** The path to `simpleGWAS_VV_Template.sh` inside the `cp` command is hard-coded. Update it to your local copy before running.

---

### `src/runFaSTLMM.py`

**Version 4.0** — Core Python script wrapping FaST-LMM's `single_snp` function. Supports multiple phenotype files, batched permutation testing, structured logging per phenotype, and automatic cleanup of FaST-LMM temporary files.

**Usage:**
```bash
python runFaSTLMM.py \
  -g <genotype_plink_prefix> \
  -p <pheno1.norm.phen> [<pheno2.norm.phen> ...] \
  -o <output_directory> \
  [-k <kinship_plink_prefix>] \
  [-c <covariance_file>] \
  [-n <nb_permutations>] \
  [-t <nb_threads>] \
  [--double-id]
```

**Arguments:**

| Flag | Description | Default |
|------|-------------|---------|
| `-g` / `--genotype` | PLINK genotype prefix | required |
| `-p` / `--phenotypes` | One or more `.norm.phen` files | required |
| `-o` / `--outdir` | Output root directory | required |
| `-k` / `--kinship` | PLINK kinship prefix (`.bed` appended automatically) | `""` |
| `-c` / `--covariance` | Covariance matrix file | `""` |
| `-n` / `--nb-permutations` | Permutation count; `0` disables permutation testing | `100` |
| `-t` / `--threads` | Threads for `LocalMultiProc` runner | `1` |
| `--double-id` | Set FID = IID = full sample name (default: split on `_`) | off |

**Key implementation details:**

- **Sample ID handling:** By default, sample names containing `_` are split into `(FID, IID)` pairs to mirror PLINK 1.9 default behaviour. Use `--double-id` if your sample names should not be split.
- **Sample overlap logging:** At startup, warns about samples present in the phenotype file but absent from the genotype `.fam` file, or vice versa.
- **Batching:** FaST-LMM is limited to 101 phenotypes per call (`BATCH_SIZE = 101`: 1 real + up to 100 permutation columns). Larger permutation counts are automatically split into multiple batches.
- **Threshold:** Computed as the 5th percentile of the best (minimum) p-value across all permutation columns, pooled over batches.
- **Output directory:** A subdirectory `<outdir>/<phen_name>/` is created for each phenotype. If it already exists, the script raises an error rather than overwriting.
- **Cleanup:** Removes FaST-LMM's `runs/` and `.work*/` temporary directories after each phenotype run.
- **Logging:** A per-phenotype `.log` file is written with timestamps, sample overlap stats, per-batch timing, threshold value, and significant SNP count.

**Outputs per phenotype** (in `<outdir>/<phen_name>/`):

| File | Description |
|------|-------------|
| `<phen_name>.first_assoc.txt.gz` | Full association results for all SNPs (gzipped TSV) |
| `<phen_name>.threshold.txt` | Permutation threshold — format: `x\n<float>` |
| `<phen_name>.signif_snps.txt` | Rows from `first_assoc` where `PValue < threshold` |
| `<phen_name>.log` | Detailed run log |

---

### `src/rank_based_Inverse_Normal_Transformation.R`

Normalises raw phenotype values using the rank-based inverse normal transformation (INT), which is the required pre-processing step before FaST-LMM (which assumes Gaussian residuals).

**Formula** (Blom 1958, offset `c = 3/8`):
```
INT(x) = qnorm( (rank(x) - 3/8) / (N - 2×(3/8) + 1) )
```
Ties are resolved by the `"average"` method. `NA` values are excluded from the rank calculation but re-inserted into the output (before the final `na.omit` removes them).

**Usage:**
```bash
Rscript rank_based_Inverse_Normal_Transformation.R <phenotype.phen>
```

**Input** — tab-separated, with header, two columns:
```
Strain      TraitName
StrainAAA   1563.2
StrainBBB   NA
StrainCCC   4218.7
```

**Output** — written to the same directory as the input, with `.phen` replaced by `.norm.phen`. `NA` rows are removed. Format: two columns (`Strain`, normalised values), tab-separated, with header, no row names.

---

### `src/calc_GIF.R`

Computes the **Genomic Inflation Factor (λ)** from FaST-LMM association results. λ is the ratio of the median observed chi-squared statistic to the theoretical chi-squared median under the null (0.4549). A value near 1.0 indicates well-calibrated p-values; values substantially above 1.0 suggest inflation, e.g. from unaccounted population structure.

**Usage:**
```bash
Rscript calc_GIF.R <phen_name>.first_assoc.txt.gz
```

**Input:** The gzipped association results from `runFaSTLMM.py`. Rows with `NA` in column 6 (the p-value column) are removed before computing χ².

**Output:** A tab-separated file `<input_file>.lgc` with two columns: the input filename and the λ value.

---

### `src/qqman_script.R`

Generates **Manhattan and QQ plots** from association results using the R `qqman` package. Supports three GWAS software output formats selected at runtime.

**Usage:**
```bash
Rscript qqman_script.R \
  <phen_name>.first_assoc.txt.gz \
  <phen_name>.threshold.txt \
  [<causal_snps_file>] \
  <gwas_soft>
```

**Arguments:**

| Position | Description |
|----------|-------------|
| `args[1]` | Association results file |
| `args[2]` | Threshold file from `runFaSTLMM.py`, or `"0"` for no threshold line |
| `args[3]` | (Optional) File of known/causal SNP IDs to highlight; if the file does not exist or is omitted, this argument is treated as `gwas_soft` |
| `args[4]` | GWAS software label: `fastlmm207`, `fastlmm`, or `gemma` |

**`fastlmm207` mode** (current pipeline default):
- Reads first 11 columns; renames to: `sid_index`, `SNP`, `CHR`, `GenDist`, `BP`, `P`, `SnpWeight`, `SnpWeightSE`, `SnpFractVarExpl`, `Mixing`, `Nullh2`.
- Manhattan plot uses an alternating two-colour scheme (`cornflowerblue` / `gray40`) with a genome-wide line at `-log10(threshold)`.
- Optionally highlights causal/known SNPs.
- Outputs: `<input>_p.jpg` (Manhattan) and `<input>_qq.jpg` (QQ plot). An R warning log is saved to `<input>_R.log`.

**`gemma` mode:**
Generates separate plots for WALD, LRT, and SCORE statistics. With a causal SNP file, also generates ROC curves (requires the `pROC` package). Outputs EPS + JPEG for each test.

**`fastlmm` mode:**
Older FaST-LMM format — single p-value column, annotated with the permutation threshold. Outputs EPS Manhattan and JPEG QQ plot.

---

### `src/addLinkageGroups.py`

Annotates significant SNPs with **linkage group** labels using a graph-based approach over PLINK LD results. Each connected component in the LD graph (where edges link SNP pairs exceeding the r² threshold within a genomic window) becomes one linkage group.

**Usage:**
```bash
python addLinkageGroups.py \
  -g <genotype_plink_prefix> \
  -r <phen_name>.signif_snps.txt \
  [-w <window_size_kb>] \
  [-r2 <r2_threshold>] \
  [-p SNP InDel SV CNV]
```

**Arguments:**

| Flag | Description | Default |
|------|-------------|---------|
| `-g` / `--genotype` | PLINK genotype prefix | required |
| `-r` / `--results` | Significant SNPs file (`*.signif_snps.txt`) | required |
| `-w` / `--window_size` | Max distance (kb) between two variants for LD consideration | `50` |
| `-r2` / `--r2` | r² threshold for placing two variants in the same group | `0.5` |
| `-p` / `--prefix_type` | Space-separated variant-type prefixes (e.g. `SNP InDel SV CNV`); triggers per-type grouping | none |

**How it works:**
1. Writes significant SNP IDs to a temporary file and calls PLINK (`--r2`, `--ld-window-r2`, `--ld-window-kb`) to compute pairwise LD.
2. Builds a `networkx` undirected graph — nodes are SNPs, edges connect pairs above the r² threshold.
3. Labels each connected component `LD_1`, `LD_2`, … and stores the result in a `LinkageGroup` column.
4. If `--prefix_type` is given, repeats graph construction within each variant type (filtering both SNPs and LD edges by prefix) and adds a `LinkageGroup_PerVariantType` column (e.g. `LD_SNP_1`, `LD_InDel_2`).
5. Removes all temporary PLINK intermediate files.

**Edge cases:**
- Exactly 1 significant SNP → assigned `LD_1` directly without invoking PLINK.
- 0 significant SNPs → `LinkageGroup` set to `None`.

**Output:** `<results_basename>.LinkageGroups.txt` — the original significant SNPs table extended with `LinkageGroup` (and optionally `LinkageGroup_PerVariantType`) columns, tab-separated.

---

## Input File Formats

### Genotypes
PLINK binary format: `.bed`, `.bim`, `.fam` sharing the same prefix. Pass only the prefix (no extension) to all scripts.

### Phenotype file (`.phen`)
Tab-separated, **with header**. Two columns: strain/sample ID and trait values. `NA` values are accepted and removed during normalisation.

```
Strain      TraitName
StrainAAA   1563.2
StrainBBB   4218.7
StrainCCC   NA
```

Place all `.phen` files for a run in the same `PHENODIR` directory.

### Normalised phenotype file (`.norm.phen`)
Same two-column format as `.phen` but with rank-based INT values and no `NA` rows. Generated automatically by the pipeline; can also be provided pre-computed to skip normalisation.

### Kinship matrix (optional)
PLINK binary format (`.bed`/`.bim`/`.fam`). If omitted, FaST-LMM constructs the GRM from the test genotypes.

### Covariance matrix (optional)
Tab-separated, no header. Columns: FID, IID, covariate value(s). Follows the PLINK alternate phenotype/covariate file convention.

---

## Output Files

For each phenotype `<COND>`, results are written to a timestamped directory `<WORKDIR>/<YYYYMMDD_HHMMSS>_fastlmm_results/<COND>/`:

| File | Description |
|------|-------------|
| `<COND>.first_assoc.txt.gz` | Full FaST-LMM association results for all tested SNPs |
| `<COND>.threshold.txt` | Permutation-based significance threshold |
| `<COND>.signif_snps.txt` | SNPs with `PValue < threshold` |
| `<COND>.signif_snps.LinkageGroups.txt` | Significant SNPs annotated with linkage group labels |
| `<COND>.first_assoc.txt.gz.lgc` | Genomic Inflation Factor (λ) |
| `<COND>.first_assoc.txt.gz_p.jpg` | Manhattan plot |
| `<COND>.first_assoc.txt.gz_qq.jpg` | QQ plot |
| `<COND>.first_assoc.txt.gz_p.eps` | Manhattan plot (vector/PostScript format) |
| `<COND>.log` | Per-phenotype run log with timing and sample overlap info |

---

## References

- Lippert et al. (2011). *FaST linear mixed models for genome-wide association studies.* Nature Methods, 8:833–835. https://doi.org/10.1038/nmeth.1681
- Blom, G. (1958). *Statistical Estimates and Transformed Beta-Variables.* Wiley.
- Beasley et al. (2009). *Rank-based inverse normal transformations are increasingly used, but are they merited?* Behavior Genetics, 39(5):580–595.