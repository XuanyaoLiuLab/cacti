
# CACTI: Leveraging Correlated Regulatory Elements for Powerful Chromatin QTL Detection

## Overview

<figure>
<img src="man/figures/cacti_overview.png" alt="CACTI overview" />
<figcaption aria-hidden="true">CACTI overview</figcaption>
</figure>

CACTI (**C**hrom**A**tin quantitative lo**C**i of mul**TI**ple peaks)
implements the *peak-window* version of the CACTI method for chromatin
QTL mapping.

The package provides a pipeline to:

- group nearby peaks into **non-overlapping windows** of fixed genomic
  size,
- residualize peak intensities with respect to covariates,
- run **multi-peak** association tests using a principal-component
  omnibus (PCO) test, and
- aggregate results across chromosomes and compute **window-level FDR**.

------------------------------------------------------------------------

## Installation

``` r
install.packages("devtools")   # if not installed yet
devtools::install_github("liliw-w/cacti", build_vignettes = TRUE)
```

After installation:

``` r
library(cacti)
```

------------------------------------------------------------------------

## Usage

After installing the package, see [full documentation and
vignettes](https://liliw-w.github.io/cacti/).

Also see in the package -

- The vignette **“CACTI peak-window pipeline”** that gives a more
  detailed walk through of inputs, outputs, and example usage.

``` r
vignette("cacti_peak_window", package = "cacti")
```

- Main functions documentation -

``` r
?cacti
?cacti_run_chr
?cacti_run_genome
?cacti_add_fdr
```

------------------------------------------------------------------------

## Bundled example data

The package includes a small toy dataset derived from Aracena et al. to
demonstrate the pipeline:

``` text
inst/extdata/
  test_pheno_meta.bed
  test_pheno.txt
  test_covariates.txt
  test_qtl_sum_stats_chr5.txt.gz
```

You can access these files using `system.file()`:

``` r
file_pheno_meta <- system.file("extdata", "test_pheno_meta.bed", package = "cacti")

file_pheno <- system.file("extdata", "test_pheno.txt", package = "cacti")

file_cov <- system.file("extdata", "test_covariates.txt", package = "cacti")

qtl_file <- system.file("extdata", "test_qtl_sum_stats_chr5.txt.gz", package = "cacti")
```

------------------------------------------------------------------------

## Input file formats

### Peak metadata (`file_pheno_meta`)

A BED-like table describing peaks:

- Example: `test_pheno_meta.bed`
- Required columns:
  - `phe_chr`: chromosome (e.g. `"chr5"`)
  - `phe_from`: peak start (integer, `< phe_to`)
  - `phe_to`: peak end (integer)
  - `phe_id`: peak identifier (unique)

### Phenotype matrix (`file_pheno`)

Peak intensity matrix:

- Example: `test_pheno.txt`
- Layout:
  - Rows = features (peaks)
  - First column = feature ID (matching `phe_id`)
  - Remaining columns = samples

### Covariate matrix (`file_cov`)

Covariates per sample:

- Example: `test_covariates.txt`
- Layout:
  - Rows = covariates
  - First column = covariate ID
  - Remaining columns = samples

### QTL summary statistics (`qtl_file`)

cis-QTL summary statistics per chromosome:

- Example: `test_qtl_sum_stats_chr5.txt.gz`
- Required columns:
  - `phe_id`: peak ID (matching phenotype and metadata)
  - `var_id`: variant ID
  - `z`: z-score for the peak–SNP association

------------------------------------------------------------------------

## Quick start: single-chromosome pipeline (`cacti_run_chr()`)

`cacti_run_chr()` runs the full peak-window pipeline for **one
chromosome**:

``` r
out_prefix <- tempfile("cacti_chr5_")

res_single_chr <- cacti_run_chr(
  window_size     = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno      = file_pheno,
  file_cov        = file_cov,
  chr             = "chr5",
  qtl_file        = qtl_file,
  out_prefix      = out_prefix,
  min_peaks       = 2
)

res_single_chr
```

This will:

1.  Group peaks into non-overlapping windows (50 kb).  
2.  Residualize the phenotype matrix with respect to covariates.  
3.  Run CACTI association tests on chr5.

`res_single_chr` is a list containing the paths to:

- `file_peak_group` – window-level group table  
- `file_peak_group_peaklevel` – peak-level group assignment  
- `file_pheno_cov_residual` – residualized phenotype matrix  
- `file_p_peak_group` – window-level p-values for chr5

------------------------------------------------------------------------

## Quick start: genome-wide wrapper (`cacti_run_genome()`)

`cacti_run_genome()` loops over chromosomes, runs `cacti_run_chr()` for
each, and then combines results to compute window-level FDR.

For the toy data, we only have chr5, but the usage pattern is the same:

``` r
out_prefix <- tempfile("cacti_genome_")

res_genome <- cacti_run_genome(
  window_size     = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno      = file_pheno,
  file_cov        = file_cov,
  chrs            = "chr5",                 # toy example: one chromosome
  qtl_files       = qtl_file,
  out_prefix      = out_prefix,
  file_fdr_out    = file.path(tempdir(), "cacti_fdr_chr5.txt.gz")
)

res_genome
```

`res_genome` is a list containing the paths to:

- `file_peak_group` – window-level group table  
- `file_peak_group_peaklevel` – peak-level group assignment  
- `file_pheno_cov_residual` – residualized phenotype matrix  
- `file_p_peak_group` – window-level p-value files for all chromosome
- `file_fdr_out` – combined FDR-added window-level results.

## Citation

If you use the CACTI method, please cite:

> Wang, L., & Liu, X. (2025). Improved chromatin QTL mapping with CACTI.
> bioRxiv, 2025-06.
