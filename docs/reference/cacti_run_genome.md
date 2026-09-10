# Run CACTI peak-window pipeline genome-wide

This function is a convenience wrapper that runs the CACTI peak-window
pipeline across multiple chromosomes and optionally computes
window-level FDR across all peak windows and chromosomes.

## Usage

``` r
cacti_run_genome(
  window_size,
  file_pheno_meta,
  file_pheno,
  file_cov,
  chrs,
  qtl_files = NULL,
  file_vcf = NULL,
  file_geno = NULL,
  file_snp_pos = NULL,
  cis_dist = 1e+05,
  p_threshold = 1,
  out_prefix,
  dir_pco = system.file("pco", package = "cacti"),
  min_peaks = 1,
  file_fdr_out = NULL,
  do_fdr = TRUE
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

- chrs:

  Character vector of chromosome labels (e.g., `paste0("chr", 1:22)`).

- qtl_files:

  Either:

  - a character vector of the same length as `chrs`, giving the cis-QTL
    file path for each chromosome in order; or

  - a single template string containing the placeholder `"{chr}"`, e.g.,
    `"extdata/test_qtl_sum_stats_{chr}.txt.gz"`. In that case, the
    placeholder is replaced by each element of `chrs`.

  If `NULL`, MatrixEQTL is run once first from genotype + phenotype +
  covariates to generate a CACTI-compatible cis-QTL file used for all
  `chrs`.

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

- file_fdr_out:

  Optional output path for the FDR-added window-level file. If `NULL`, a
  default filename is constructed from `out_prefix` and `window_size`.

- do_fdr:

  Logical; if `TRUE` (default), run
  [`cacti_add_fdr()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_add_fdr.md)
  across all chromosomes. If `FALSE`, skip FDR correction.

## Value

Invisibly returns a named list of output paths with elements:

- file_peak_group:

  Path to the window-level group.

- file_peak_group_peaklevel:

  Path to the peak-level group.

- file_pheno_residual:

  Path to the residualized phenotype matrix.

- file_p_peak_group:

  Path to the per-window p-value file for all chromosome.

- file_fdr_out:

  Path to the FDR-added window-level result file (`NULL` when
  `do_fdr = FALSE`).

## Details

For each chromosome in `chrs`, it runs per-window p-value calculation
and collects per-chromosome p-value files. If `do_fdr = TRUE`, it then
calls
[`cacti_add_fdr()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_add_fdr.md)
once, aggregating all chromosomes to obtain q-values for the top-hit
p-values in each window.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example (chr5 as all chr)

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

out_prefix <- tempfile("cacti_genome_")

res <- cacti_run_genome(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  chrs = "chr5",
  qtl_files = qtl_file,
  out_prefix = out_prefix,
  dir_pco = system.file("pco", package = "cacti"),
  min_peaks = 1,
  file_fdr_out = file.path(tempdir(), "cacti_fdr_chr5.txt.gz")
)
} # }
```
