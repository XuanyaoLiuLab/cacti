# Run CACTI peak-window pipeline end-to-end

High-level wrapper for the CACTI peak-window workflow. This function
supports two input modes: 1) genotype + phenotype + covariates (runs
MatrixEQTL first), 2) precomputed summary statistics (skips MatrixEQTL).

## Usage

``` r
cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta,
  file_pheno,
  file_cov,
  out_prefix,
  chr = "All",
  qtl_files = NULL,
  qtl_file = NULL,
  file_vcf = NULL,
  file_geno = NULL,
  file_snp_pos = NULL,
  cis_dist = 1e+05,
  dir_pco = system.file("pco", package = "cacti"),
  min_peaks = 1,
  do_fdr = TRUE,
  file_fdr_out = NULL
)
```

## Arguments

- window_size:

  Size of non-overlapping genomic windows used to group peaks before
  CACTI testing. Accepts values like `"50kb"` (kilobases) or numeric
  base-pair values (e.g., `50000`).

- file_pheno_meta:

  Path to phenotype meta table with columns: `phe_id`, `phe_chr`,
  `phe_from`, `phe_to`.

- file_pheno:

  Path to phenotype matrix. Rows = features, first column = ID.

- file_cov:

  Path to covariate matrix. Rows = covariates, first column = ID.

- out_prefix:

  Output prefix.

- chr:

  Chromosome selection. Use `"All"` (default) to iterate over all
  chromosomes in `file_pheno_meta`, or provide one chromosome label
  (e.g., `"chr5"`), or a character vector of labels.

- qtl_files:

  Optional summary-stat input (vector, single path, or `{chr}`
  template). If `NULL`, MatrixEQTL is run first from genotype input.

- qtl_file:

  Optional single-chromosome alias of `qtl_files`.

- file_vcf:

  Optional VCF genotype input for genotype-input mode.

- file_geno:

  Optional genotype matrix input for genotype-input mode (if no VCF).

- file_snp_pos:

  Optional SNP-position input for genotype-input mode (if no VCF).

- cis_dist:

  Cis distance for MatrixEQTL (default 100kb).

- dir_pco:

  Path to PCO helper files.

- min_peaks:

  Minimum number of peaks required for a window to be included in
  testing. Included windows with 1 peak use univariate p-values;
  included windows with \>=2 peaks use PCO.

- do_fdr:

  Logical; if `TRUE` (default), run FDR correction across all processed
  chromosomes. If `FALSE`, skip FDR correction.

- file_fdr_out:

  Optional output path for FDR-added result file.

## Value

Invisibly returns the output list from
[`cacti_run_genome()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_genome.md).

## Details

By default (`chr = "All"`), it runs the genome-wide CACTI peak-window
workflow via
[`cacti_run_genome()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_genome.md)
across all chromosomes present in `file_pheno_meta` and performs FDR
correction.
