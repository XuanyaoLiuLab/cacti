# Run CACTI peak-window pipeline for a single chromosome

Steps:

1.  Group peaks into non-overlapping windows.

2.  Residualize phenotype by covariates.

3.  Run association test on peak groups for a single chromosome.

## Usage

``` r
cacti_run_chr(
  window_size,
  file_pheno_meta,
  file_pheno,
  file_cov,
  chr,
  qtl_file = NULL,
  file_vcf = NULL,
  file_geno = NULL,
  file_snp_pos = NULL,
  cis_dist = 1e+05,
  p_threshold = 1,
  out_prefix,
  dir_pco = system.file("pco", package = "cacti"),
  min_peaks = 1
)
```

## Arguments

- window_size:

  Character like `"50kb"` (kb units) or numeric (bp).

- file_pheno_meta:

  Path to input phenotype meta table (with header). Columns are -

  - `phe_chr`: chromosome, with "chr" prefix (e.g. "chr1")

  - `phe_from`: peak start (integer, \< `phe_to`)

  - `phe_to`: peak end (integer)

  - `phe_id`: peak identifier (unique per row)

- file_pheno:

  Path to phenotype expression matrix. Rows = features; first column =
  feature ID; remaining columns = samples.

- file_cov:

  Path to covariate matrix. Rows = covariates; first column = covariate
  ID; remaining columns = samples.

- chr:

  Chromosome label to process (e.g., "chr1", "chr21").

- qtl_file:

  Optional path to cis QTL summary stats file for this chromosome. Must
  contain columns: phe_id, var_id, z. If `NULL`, MatrixEQTL is run first
  from genotype + phenotype + covariates to generate this file.

- file_vcf:

  Optional path to input VCF file for MatrixEQTL preprocessing.

- file_geno:

  Optional path to genotype matrix if no VCF is provided.

- file_snp_pos:

  Optional path to SNP-position file if no VCF is provided.

- cis_dist:

  Cis-window distance for MatrixEQTL (default 100000 bp).

- p_threshold:

  P-value threshold for MatrixEQTL output (default 1.0).

- out_prefix:

  Output prefix used to construct all output filenames.

- dir_pco:

  Directory containing association test helpers:
  ModifiedPCOMerged_acat.R, liu.R, liumod.R, davies.R, qfc.so.

- min_peaks:

  Minimum number of peaks required for a window to be included in
  testing. Included windows with 1 peak use univariate p-values;
  included windows with \>=2 peaks use PCO.

## Value

Invisibly returns a named list of output paths with elements:

- file_peak_group:

  Path to the window-level group.

- file_peak_group_peaklevel:

  Path to the peak-level group.

- file_pheno_residual:

  Path to the residualized phenotype matrix.

- file_p_peak_group:

  Path to the per-window p-value file for this chromosome.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example (chr5 only)

file_pheno_meta <- system.file(
  "extdata", "test_pheno_meta.bed",
  package = "cacti"
)

file_pheno <- system.file(
  "extdata", "test_pheno.txt",
  package = "cacti"
)

file_cov <- system.file(
  "extdata", "test_covariates.txt",
  package = "cacti"
)

qtl_file <- system.file(
  "extdata", "test_qtl_sum_stats_chr5.txt.gz",
  package = "cacti"
)

out_prefix <- tempfile("cacti_chr5_")

res <- cacti_run_chr(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  chr = "chr5",
  qtl_file = qtl_file,
  out_prefix = out_prefix,
  dir_pco = system.file("pco", package = "cacti"),
  min_peaks = 1
)
} # }
```
