# CACTI peak-window pipeline

## Introduction

This vignette explains:

1.  The expected **input file formats**.  
2.  Running CACTI.  
3.  The **output files** produced by the pipeline.

For those who want an in-depth understanding of the workflow, we explain:

1.  The main functions.
2.  The overall workflow.  
3.  Example usage with bundled toy data under `inst/extdata/`.

Two input modes are available:

- **Mode A: genotype + phenotype + covariates**  
  CACTI takes standard genotype, phenotype, and covariate input, and automatically maps QTLs of multiple peaks in pre-defined windows.
- **Mode B: precomputed summary statistics**  
  If summary statistics of univariate associations (each SNP and each peak) are available, CACTI takes the summary statistics and performs window-based QTL calling.

By default, `cacti_peak_window()` runs **genome-wide** (`chr = "All"`),
iterating through all chromosomes present in `file_pheno_meta` and performing
FDR correction.

------------------------------------------------------------------------

## 1) Expected input file formats

### Required inputs (mode A: genotype + phenotype + covariates)

1.  **Genotype**

- `file_vcf`, or  
- `file_geno` + `file_snp_pos`

2.  **Phenotype matrix** (`file_pheno`)

- rows = peaks/features  
- first column = feature ID (`ID`)  
- remaining columns = sample IDs  
- should include features from all chromosomes for genome-wide runs

3.  **Phenotype metadata** (`file_pheno_meta`)

- required columns:
    - `phe_id`
    - `phe_chr`
    - `phe_from`
    - `phe_to`  
- should include all chromosomes to be analyzed

4.  **Covariate matrix** (`file_cov`)

- rows = covariates  
- first column = covariate ID  
- remaining columns = sample IDs

For genome-wide analysis, keep genotype, phenotype, and metadata as
all-chromosome inputs. Chromosome iteration is controlled by `chr`.

These inputs follow standard QTL mapping conventions used by MatrixEQTL and
QTLtools (genotype, phenotype, and covariates with matched sample IDs), with
CACTI-specific phenotype metadata columns (`phe_id`, `phe_chr`, `phe_from`,
`phe_to`) used for window-based grouping and chromosome iteration.

``` r
library(cacti)

file_pheno_meta <- system.file("extdata", "test_cacti_peak_chr5_pheno_meta.bed", package = "cacti")
file_pheno <- system.file("extdata", "test_cacti_peak_chr5_pheno.txt", package = "cacti")
file_cov <- system.file("extdata", "test_cacti_peak_chr5_covariates.txt", package = "cacti")
file_vcf <- system.file("extdata", "test_cacti_peak_chr5_geno.vcf", package = "cacti")

head(data.table::fread(file_pheno_meta))
head(data.table::fread(file_pheno))[, 1:5]
head(data.table::fread(file_cov))[, 1:5]
```

### Optional input (summary-statistics mode)

If you already have cis-QTL summary stats, you can provide:

- `qtl_files` (or `qtl_file`) with required columns:
    - `phe_id`
    - `var_id`
    - `z`

``` r
qtl_file <- system.file("extdata", "test_cacti_peak_chr5_matrixqtl_sumstats.txt.gz", package = "cacti")
head(data.table::fread(qtl_file))
```

------------------------------------------------------------------------

## 2) Running CACTI

Use one function:
[`cacti_peak_window()`](https://liliw-w.github.io/cacti/reference/cacti_peak_window.md).

`window_size` controls how peaks are grouped into non-overlapping genomic windows
before testing (for example, `"50kb"` means 50,000-bp windows).

### A. Mode A: genotype + phenotype + covariates

This runs MatrixEQTL first, then standard CACTI automatically.

``` r
library(cacti)

res <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  file_vcf = file_vcf,
  cis_dist = 100000,         # 100 kb cis window for MatrixEQTL
  chr = "All",               # default; can be omitted
  do_fdr = TRUE,             # default; can be omitted
  out_prefix = tempfile("cacti_peak_window_")
)

res
```

### B. Mode B: summary stats already available

If you already have QTL summary statistics, pass `qtl_files` (or `qtl_file`).
In this mode, MatrixEQTL is skipped.

``` r
res_summary <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  chr = "All",
  do_fdr = TRUE,
  qtl_file = qtl_file,       # optional summary-stat mode
  out_prefix = tempfile("cacti_peak_window_summary_")
)

res_summary
```

------------------------------------------------------------------------

## 3) Output files produced by the pipeline

`cacti_peak_window()` returns the same output structure as
[`cacti_run_genome()`](https://liliw-w.github.io/cacti/reference/cacti_run_genome.md):

- `file_peak_group`
  - grouped windows table
- `file_peak_group_peaklevel`
  - peak-to-window mapping
- `file_pheno_cov_residual`
  - residualized phenotype matrix
- `file_p_peak_group`
  - per-chromosome `(group, snp, pval)` files
- `file_fdr_out`
  - final FDR-added window-level results

``` r
names(res_summary)
res_summary$file_fdr_out
```

For detailed low-level APIs, see:

- [`?cacti_peak_window`](https://liliw-w.github.io/cacti/reference/cacti_peak_window.md)  
- [`?cacti_run_chr`](https://liliw-w.github.io/cacti/reference/cacti_run_chr.md)  
- [`?cacti_run_genome`](https://liliw-w.github.io/cacti/reference/cacti_run_genome.md)  
- [`?cacti_matrixqtl_cis`](https://liliw-w.github.io/cacti/reference/cacti_matrixqtl_cis.md)

------------------------------------------------------------------------

## 4) In-depth: Main functions

### High-level entry

- `cacti_peak_window()`
  - recommended one-function entry point.  
  - mode A: genotype input (`file_vcf` or `file_geno` + `file_snp_pos`) to run MatrixEQTL first, then CACTI.  
  - mode B: summary-stat input (`qtl_file` / `qtl_files`) to skip MatrixEQTL.  
  - default chromosome behavior: `chr = "All"` (iterate over all chromosomes in `file_pheno_meta`).  
  - optional `do_fdr` controls FDR correction.

### Core wrappers

- `cacti_run_chr()`
  - optional low-level wrapper for targeted single-chromosome runs.  
  - can generate QTL stats first when `qtl_file = NULL`.

- `cacti_run_genome()`
  - runs CACTI across multiple chromosomes and adds FDR.  
  - can generate QTL stats first when `qtl_files = NULL`.

### Supporting steps

- `cacti_matrixqtl_cis()`: runs MatrixEQTL and writes CACTI-compatible summary stats.
- `cacti_group_peak_window()`: groups peaks into windows.
- `cacti_pheno_cov_residual()`: residualizes phenotypes by covariates.
- `cacti_cal_p()`: computes per-window p-values.
- `cacti_add_fdr()`: aggregates p-values and adds FDR.

------------------------------------------------------------------------

## 5) In-depth: Overall workflow

### Mode A: genotype + phenotype + covariates

1.  Input all-chromosome genotype (`file_vcf` or `file_geno` + `file_snp_pos`), phenotype, metadata, and covariates.  
2.  Run MatrixEQTL to produce cis-QTL summary stats (`phe_id`, `var_id`, `z`, `pval`).  
3.  Group peaks into non-overlapping windows.  
4.  Residualize phenotype by covariates.  
5.  For each window/SNP:
    - use univariate test for 1-peak windows,
    - use PCO test for windows with >=2 peaks.
6.  By default (`chr = "All"`), iterate over all chromosomes and aggregate with FDR.

### Mode B: precomputed summary statistics

1.  Provide precomputed summary stats (`qtl_file`/`qtl_files`) for all chromosomes you want to run.  
2.  Skip MatrixEQTL and directly run CACTI steps 3-6 above.

------------------------------------------------------------------------

## 6) In-depth: Example usage with bundled toy data

### A. One-function workflow (recommended)

``` r
library(cacti)

file_pheno_meta <- system.file("extdata", "test_cacti_peak_chr5_pheno_meta.bed", package = "cacti")
file_pheno <- system.file("extdata", "test_cacti_peak_chr5_pheno.txt", package = "cacti")
file_cov <- system.file("extdata", "test_cacti_peak_chr5_covariates.txt", package = "cacti")
file_vcf <- system.file("extdata", "test_cacti_peak_chr5_geno.vcf", package = "cacti")

res <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  file_vcf = file_vcf,
  chr = "All",
  do_fdr = TRUE,
  cis_dist = 100000,
  out_prefix = tempfile("cacti_peak_window_in_depth_")
)

res
```

### B. Optional summary-stat workflow

``` r
qtl_file <- system.file("extdata", "test_cacti_peak_chr5_matrixqtl_sumstats.txt.gz", package = "cacti")

res_qtl <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  qtl_file = qtl_file,
  chr = "All",
  do_fdr = TRUE,
  out_prefix = tempfile("cacti_peak_window_in_depth_qtl_")
)

res_qtl
```

### C. Low-level wrapper workflow

``` r
# Genome-wide wrapper (+ FDR)
res_genome <- cacti_run_genome(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  chrs = paste0("chr", 1:22),
  qtl_files = qtl_file,
  out_prefix = tempfile("cacti_run_genome_"),
  min_peaks = 1,
  do_fdr = TRUE
)

# Optional: single chromosome in the same high-level API
res_chr <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  file_vcf = file_vcf,
  chr = "chr5",
  do_fdr = FALSE,  # or TRUE
  out_prefix = tempfile("cacti_peak_window_chr5_")
)
```
