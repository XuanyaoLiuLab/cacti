# Run Nominal Cis-QTL Mapping (CACTI-S)

Performs linear regression using MatrixEQTL to generate summary
statistics.

## Usage

``` r
cacti_s_map_cis(
  file_pheno,
  file_pheno_meta,
  file_cov,
  file_qtl_out,
  file_vcf = NULL,
  file_geno = NULL,
  file_snp_pos = NULL,
  cis_dist = 1e+05,
  p_threshold = 1
)
```

## Arguments

- file_pheno:

  Path to the processed phenotype file with position meta and expression
  data (from `cacti_s_preprocess`).

- file_pheno_meta:

  Path for the meta file of processed phenotype.

- file_cov:

  Path to covariate matrix.

- file_qtl_out:

  Output path for summary statistics.

- file_vcf:

  Path to input VCF file.

- file_geno:

  (Optional) Path to genotype matrix if no VCF.

- file_snp_pos:

  (Optional) Path to SNP positions if no VCF.

- cis_dist:

  Cis-window distance (default 100000 bp = 100kb).

- p_threshold:

  P-value threshold for output (default 1.0, print all associations).

## Value

Invisibly returns the summary stats.

## Examples

``` r
if (FALSE) { # \dontrun{
# Use local paths for testing
file_pheno <- "inst/extdata/test_results/test_pheno_norm.txt"
file_pheno_meta <- "inst/extdata/test_results/test_pheno_norm_meta.txt"
file_cov   <- "inst/extdata/test_cov.txt"
file_vcf   <- "inst/extdata/test_geno.vcf"

file_out <- "inst/extdata/test_results/test_cis_qtl_stats.txt"

# Run
if (file.exists(file_pheno) && file.exists(file_vcf)) {
  cacti_s_map_cis(
    file_pheno   = file_pheno,
    file_pheno_meta = file_pheno_meta,
    file_cov     = file_cov,
    file_vcf     = file_vcf,
    file_qtl_out = file_out,
    cis_dist     = 100000
  )
}
} # }
```
