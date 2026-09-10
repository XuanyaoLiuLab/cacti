# Run MatrixEQTL and write CACTI-compatible cis-QTL summary stats

Thin wrapper around
[`cacti_s_map_cis()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_map_cis.md)
that generates a summary-statistics table with columns `phe_id`,
`var_id`, `z`, `pval` for downstream CACTI peak-window testing. The `z`
column contains signed Z-scores.

## Usage

``` r
cacti_matrixqtl_cis(
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

Invisibly returns the summary-statistics data frame.
