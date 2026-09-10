# CACTI: peak-based data preparation & mapping

## Introduction

CACTI groups nearby peaks into windows and tests their joint association
with genetic variants. This vignette describes the inputs, a
one-function quick start, the output files, and the individual analysis
steps.

Two input modes are available:

- **Genotype input:** provide genotype, phenotype, phenotype metadata,
  and covariates. CACTI runs single-peak mapping with MatrixEQTL before
  testing windows.
- **Summary-statistics input:** provide precomputed SNP–peak association
  statistics instead of genotype. Phenotypes and covariates are still
  required to estimate the correlation between peaks.

The main function is
[`cacti_peak_window()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_peak_window.md).
By default, `chr = "All"` processes all chromosomes present in the
phenotype metadata, and `do_fdr = TRUE` adds window-level q-values.

## Part 1: Overview of the CACTI pipeline

### Workflow

1.  **Map single-peak cis-QTLs.** Run MatrixEQTL on genotype, phenotype,
    and covariates. Skip this step when summary statistics are supplied.
2.  **Group peaks into windows.** Assign peaks to non-overlapping
    genomic windows using their start positions.
3.  **Adjust phenotypes for covariates.** Use the residualized
    phenotypes to estimate correlations between peaks.
4.  **Test windows and add FDR.** Windows with two or more peaks use the
    principal-component omnibus test. Included single-peak windows use a
    univariate test. Aggregate the results and add window-level
    q-values.

| Function | Role |
|:---|:---|
| [`cacti_peak_window()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_peak_window.md) | Runs the complete workflow. |
| [`cacti_matrixqtl_cis()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_matrixqtl_cis.md) | Generates single-peak cis-QTL summary statistics. |
| [`cacti_group_peak_window()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_group_peak_window.md) | Groups peaks into windows. |
| [`cacti_pheno_cov_residual()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_pheno_cov_residual.md) | Adjusts phenotypes for covariates. |
| [`cacti_cal_p()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_cal_p.md) | Tests SNP–window associations for one chromosome. |
| [`cacti_add_fdr()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_add_fdr.md) | Summarizes window results and adds q-values. |

### Input files

| Input | Format |
|:---|:---|
| `file_pheno_meta` | Four columns, in this order: `phe_id`, `phe_chr`, `phe_from`, `phe_to`. Use consistent chromosome labels, such as `chr5`. |
| `file_pheno` | Rows are peaks; the first column is the peak ID and the remaining columns are sample IDs. Supply normalized peak intensities. |
| `file_cov` | Rows are covariates; the first column is the covariate ID and the remaining columns are sample IDs. |
| `file_vcf` | Genotypes in VCF format. Alternatively, supply `file_geno` and `file_snp_pos`; see [`?cacti_matrixqtl_cis`](https://xuanyaoliulab.github.io/cacti/reference/cacti_matrixqtl_cis.md). |
| `qtl_file` or `qtl_files` | Alternative to genotype input: a table containing `phe_id`, `var_id`, and `z` (signed Z-scores). Include the associations to be tested, not only significant hits. |

Use matching peak and sample IDs across inputs. In summary-statistics
mode, the association statistics, phenotypes, and covariates should
describe the same analysis samples and features. For genome-wide
analysis, provide the inputs for all chromosomes to be analyzed.

## Part 2: Quick start

### 1. Locate test data

The package includes a small chromosome 5 example. All inputs below are
read from the installed package; results are written to a temporary
directory.

``` r

library(cacti)

file_pheno_meta <- system.file("extdata", "test_cacti_peak_chr5_pheno_meta.bed", package = "cacti")
file_pheno <- system.file("extdata", "test_cacti_peak_chr5_pheno.txt", package = "cacti")
file_cov <- system.file("extdata", "test_cacti_peak_chr5_covariates.txt", package = "cacti")
file_vcf <- system.file("extdata", "test_cacti_peak_chr5_geno.vcf", package = "cacti")
qtl_file <- system.file("extdata", "test_cacti_peak_chr5_matrixqtl_sumstats.txt.gz", package = "cacti")

stopifnot(all(file.exists(c(file_pheno_meta, file_pheno, file_cov, file_vcf, qtl_file))))
print(basename(c(file_pheno_meta, file_pheno, file_cov, file_vcf, qtl_file)))
#> [1] "test_cacti_peak_chr5_pheno_meta.bed"           
#> [2] "test_cacti_peak_chr5_pheno.txt"                
#> [3] "test_cacti_peak_chr5_covariates.txt"           
#> [4] "test_cacti_peak_chr5_geno.vcf"                 
#> [5] "test_cacti_peak_chr5_matrixqtl_sumstats.txt.gz"

out_dir <- tempfile("cacti_peak_vignette_")
dir.create(out_dir)
out_prefix <- file.path(out_dir, "quickstart")
```

#### Phenotype metadata

``` r

knitr::kable(
  data.table::fread(file_pheno_meta, nrows = 5),
  caption = "Peak coordinates: first five rows."
)
```

| phe_id             | phe_chr | phe_from | phe_to |
|:-------------------|:--------|---------:|-------:|
| chr5_10036_11758   | chr5    |    10036 |  11758 |
| chr5_188815_189099 | chr5    |   188815 | 189099 |
| chr5_189860_192019 | chr5    |   189860 | 192019 |
| chr5_209144_209685 | chr5    |   209144 | 209685 |
| chr5_214353_214634 | chr5    |   214353 | 214634 |

Peak coordinates: first five rows. {.table}

#### Phenotype and covariate matrices

Only the first four samples are shown.

``` r

pheno_preview <- data.table::fread(file_pheno, nrows = 5)
cov_preview <- data.table::fread(file_cov, nrows = 5)

knitr::kable(
  pheno_preview[, 1:5, with = FALSE],
  digits = 3,
  caption = "Normalized peak intensities: first five peaks and four samples."
)
```

| ID                 | Sample1 | Sample2 | Sample3 | Sample4 |
|:-------------------|--------:|--------:|--------:|--------:|
| chr5_10036_11758   |   0.375 |   1.206 |   1.422 |   0.348 |
| chr5_188815_189099 |  -3.561 |  -0.159 |   0.269 |   1.417 |
| chr5_189860_192019 |   0.942 |   1.487 |  -1.165 |   1.900 |
| chr5_209144_209685 |   1.596 |  -1.086 |  -1.532 |  -1.251 |
| chr5_214353_214634 |  -0.183 |   0.600 |   0.844 |   0.908 |

Normalized peak intensities: first five peaks and four samples. {.table}

``` r

knitr::kable(
  cov_preview[, 1:5, with = FALSE],
  digits = 3,
  caption = "Covariates: first five rows and four samples."
)
```

| ID  |  Sample1 |  Sample2 | Sample3 |  Sample4 |
|:----|---------:|---------:|--------:|---------:|
| PC1 | -123.808 |   82.096 |  39.246 |   41.961 |
| PC2 |  -59.295 | -127.288 | -77.525 | -107.530 |
| PC3 |   30.379 | -106.638 | -44.708 |   40.269 |
| PC4 |  -33.170 |    3.606 | -92.398 | -126.294 |
| PC5 |    2.860 |   27.322 |  84.919 |  -72.070 |

Covariates: first five rows and four samples. {.table}

### 2. Run the pipeline

Use 50-kb windows for peak grouping and a 100-kb cis distance for
single-peak mapping. These are separate settings: `window_size` controls
peak grouping, whereas `cis_dist` controls which SNP–peak pairs are
tested by MatrixEQTL.

Identical warnings are displayed only once on this page. The example can
produce an ACAT warning when a component P value equals 1; the final
output P values are checked to be finite and within \[0, 1\].

``` r

res <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  file_vcf = file_vcf,
  cis_dist = 100000,
  chr = "All",
  min_peaks = 1,
  do_fdr = TRUE,
  out_prefix = out_prefix
)
#> Warning in ACAT::ACAT(rbind(p.PCMinP, p.PCFisher, p.PCLC, p.WI, p.Wald, : There
#> are p-values that are exactly 1!
```

In this example, `"All"` means chromosome 5 because only chromosome 5 is
present in the input metadata.

### 3. Explore results

The function returns paths to the generated files.

``` r

data.frame(
  output = names(res),
  file = vapply(res, function(x) paste(basename(x), collapse = "; "), character(1))
)
#>                                              output
#> file_peak_group                     file_peak_group
#> file_peak_group_peaklevel file_peak_group_peaklevel
#> file_pheno_cov_residual     file_pheno_cov_residual
#> file_p_peak_group                 file_p_peak_group
#> file_fdr_out                           file_fdr_out
#>                                                                       file
#> file_peak_group                       quickstart_peak_group_window50kb.txt
#> file_peak_group_peaklevel quickstart_peak_group_window50kb_peak_as_row.txt
#> file_pheno_cov_residual                  quickstart_pheno_cov_residual.txt
#> file_p_peak_group                   quickstart_pval_window50kb_chr5.txt.gz
#> file_fdr_out                   quickstart_pval_window50kb_fdr_added.txt.gz
```

| Output | Contents |
|:---|:---|
| `file_peak_group` | Window coordinates and the peaks included in each window. |
| `file_peak_group_peaklevel` | Peak-to-window mapping. |
| `file_pheno_cov_residual` | Phenotype matrix after covariate adjustment. |
| `file_p_peak_group` | SNP–window association statistics, with columns `group`, `snp`, and `pval`. |
| `file_fdr_out` | Window-level summaries and q-values. |

All previews below are read from the files generated by the quick start.

#### Peak-to-window mapping

``` r

knitr::kable(
  data.table::fread(res$file_peak_group_peaklevel, nrows = 5),
  caption = "Peak-to-window mapping: first five rows."
)
```

| phe_id             | phe_chr | phe_from | phe_to | group | n_pid_group |
|:-------------------|:--------|---------:|-------:|:------|------------:|
| chr5_10036_11758   | chr5    |    10036 |  11758 | G1    |           1 |
| chr5_188815_189099 | chr5    |   188815 | 189099 | G2    |           2 |
| chr5_189860_192019 | chr5    |   189860 | 192019 | G2    |           2 |
| chr5_209144_209685 | chr5    |   209144 | 209685 | G3    |           4 |
| chr5_214353_214634 | chr5    |   214353 | 214634 | G3    |           4 |

Peak-to-window mapping: first five rows. {.table}

#### Covariate-adjusted phenotype matrix

``` r

residual_preview <- data.table::fread(res$file_pheno_cov_residual, nrows = 5)
knitr::kable(
  residual_preview[, 1:5, with = FALSE],
  digits = 3,
  caption = "Covariate-adjusted phenotypes: first five peaks and four samples."
)
```

| ID                 | Sample1 | Sample2 | Sample3 | Sample4 |
|:-------------------|--------:|--------:|--------:|--------:|
| chr5_10036_11758   |  -0.796 |   0.102 |   0.475 |  -0.445 |
| chr5_188815_189099 |  -2.061 |  -0.548 |   0.314 |   0.885 |
| chr5_189860_192019 |   0.499 |   0.379 |  -1.227 |   0.221 |
| chr5_209144_209685 |   1.800 |  -0.912 |  -1.484 |  -1.107 |
| chr5_214353_214634 |   0.719 |   0.262 |  -0.048 |   0.461 |

Covariate-adjusted phenotypes: first five peaks and four samples.
{.table}

#### SNP–window association statistics

``` r

association_preview <- data.table::fread(res$file_p_peak_group[1], nrows = 5)
association_preview$pval <- format(association_preview$pval, digits = 3, scientific = TRUE)
knitr::kable(
  association_preview,
  caption = "SNP–window associations: first five rows, not ranked by significance."
)
```

| group | snp          | pval     |
|:------|:-------------|:---------|
| G2    | 5:180531:C:T | 9.67e-04 |
| G2    | 5:181839:A:C | 9.67e-04 |
| G2    | 5:210785:C:T | 3.41e-04 |
| G2    | 5:206988:C:A | 5.77e-04 |
| G2    | 5:207082:C:T | 5.77e-04 |

SNP–window associations: first five rows, not ranked by significance.
{.table}

#### Window-level FDR results

[`cacti_add_fdr()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_add_fdr.md)
combines SNP–window P values within each window using ACAT and
calculates q-values across the analyzed windows. The output also
includes the minimum P value, its SNP ID, the number of tested SNPs, and
a Bonferroni-adjusted minimum P value. The `q` column corresponds to
`best_hit_acat`.

``` r

fdr_preview <- data.table::fread(res$file_fdr_out, nrows = 5)
for (column in c("best_hit", "best_hit_bonf", "best_hit_acat", "q")) {
  fdr_preview[[column]] <- format(fdr_preview[[column]], digits = 3, scientific = TRUE)
}
knitr::kable(
  fdr_preview,
  caption = "Window-level summaries: first five rows."
)
```

| group | best_hit | n_var_in_cis | best_hit_bonf | best_hit_acat | snp          | q        |
|:------|:---------|-------------:|:--------------|:--------------|:-------------|:---------|
| G1    | 1.12e-02 |          563 | 1             | 9.92e-01      | 5:60134:T:A  | 9.98e-01 |
| G10   | 1.21e-02 |         1125 | 1             | 9.95e-01      | 5:632227:G:C | 9.98e-01 |
| G11   | 8.73e-04 |         1531 | 1             | 9.95e-01      | 5:945144:C:T | 9.98e-01 |
| G12   | 1.10e-02 |         1411 | 1             | 9.95e-01      | 5:853507:A:G | 9.98e-01 |
| G13   | 1.37e-02 |         1226 | 1             | 9.97e-01      | 5:975668:A:C | 9.98e-01 |

Window-level summaries: first five rows. {.table}

## Part 3: Step-by-step workflow

The calls below show the individual steps performed by the quick start.
They are shown for users who want to run or modify one step separately.
This vignette reuses the files generated above instead of rerunning the
same analysis.

### Step 1: Map single-peak cis-QTLs

The genotype-input workflow writes this intermediate summary-statistics
file:

``` r

file_qtl_generated <- paste0(out_prefix, "_matrixqtl_cis_all_chrs.txt.gz")
stopifnot(file.exists(file_qtl_generated))
```

The corresponding single-peak mapping call is:

``` r

cacti_matrixqtl_cis(
  file_pheno = file_pheno,
  file_pheno_meta = file_pheno_meta,
  file_cov = file_cov,
  file_vcf = file_vcf,
  file_qtl_out = file_qtl_generated,
  cis_dist = 100000,
  p_threshold = 1.0
)
```

`p_threshold = 1.0` retains all tested associations. The output includes
signed Z-scores (`z`) and P values (`pval`).

``` r

generated_qtl_preview <- data.table::fread(file_qtl_generated, nrows = 5)
generated_qtl_preview$pval <- format(generated_qtl_preview$pval, digits = 3, scientific = TRUE)
knitr::kable(
  generated_qtl_preview,
  digits = 4,
  caption = "Generated single-peak summary statistics: first five rows."
)
```

| phe_id               | var_id        |       z | pval     |
|:---------------------|:--------------|--------:|:---------|
| chr5_5364492_5365547 | 5:5365177:T:C |  4.4649 | 8.01e-06 |
| chr5_5116511_5117722 | 5:5119822:A:C | -4.3415 | 1.42e-05 |
| chr5_5364492_5365547 | 5:5350343:A:G |  4.3320 | 1.48e-05 |
| chr5_6221823_6222165 | 5:6318715:C:A |  4.3126 | 1.61e-05 |
| chr5_5364492_5365547 | 5:5350773:A:G |  4.2845 | 1.83e-05 |

Generated single-peak summary statistics: first five rows. {.table}

### Step 2: Group peaks into windows

``` r

cacti_group_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_peak_group = res$file_peak_group,
  file_peak_group_peaklevel = res$file_peak_group_peaklevel
)
```

The window table records the number of peaks in `n_pid_group` and their
IDs in `pid_group`. The peak-level mapping is shown above.

### Step 3: Adjust phenotypes for covariates

``` r

cacti_pheno_cov_residual(
  file_pheno = file_pheno,
  file_cov = file_cov,
  file_pheno_cov_residual = res$file_pheno_cov_residual
)
```

CACTI uses the correlations between these residualized phenotypes in the
multivariate test.

### Step 4: Test SNP–window associations

``` r

cacti_cal_p(
  file_qtl_cis_norm = file_qtl_generated,
  chr = "chr5",
  file_peak_group = res$file_peak_group,
  file_pheno_cov_residual = res$file_pheno_cov_residual,
  file_p_peak_group = res$file_p_peak_group[1],
  min_peaks = 1
)
```

`min_peaks = 1` includes both single-peak and multi-peak windows. Set it
to `2` to restrict testing to multi-peak windows. Within a multi-peak
window, the test uses SNPs with association statistics available for
every included peak.

### Step 5: Add window-level FDR

``` r

cacti_add_fdr(
  file_all_pval = res$file_p_peak_group,
  file_fdr_out = res$file_fdr_out
)
```

For a genome-wide analysis, pass all chromosome-specific P-value files
together so that q-values are calculated across the analyzed
chromosomes.

## Part 4: Optional usage

### Use precomputed summary statistics

Provide `qtl_file` instead of `file_vcf`. This skips MatrixEQTL, but the
phenotype and covariate inputs are still used for peak correlations. The
following example runs this mode using the bundled summary-statistics
file.

``` r

qtl_preview <- data.table::fread(qtl_file, nrows = 5)
qtl_preview$pval <- format(qtl_preview$pval, digits = 3, scientific = TRUE)
knitr::kable(
  qtl_preview,
  digits = 4,
  caption = "Bundled single-peak summary statistics: first five rows."
)
```

| phe_id               | var_id        |       z | pval     |
|:---------------------|:--------------|--------:|:---------|
| chr5_5364492_5365547 | 5:5365177:T:C |  4.4649 | 8.01e-06 |
| chr5_5116511_5117722 | 5:5119822:A:C | -4.3415 | 1.42e-05 |
| chr5_5364492_5365547 | 5:5350343:A:G |  4.3320 | 1.48e-05 |
| chr5_6221823_6222165 | 5:6318715:C:A |  4.3126 | 1.61e-05 |
| chr5_5364492_5365547 | 5:5350773:A:G |  4.2845 | 1.83e-05 |

Bundled single-peak summary statistics: first five rows. {.table}

``` r

res_summary <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  qtl_file = qtl_file,
  chr = "All",
  min_peaks = 1,
  do_fdr = TRUE,
  out_prefix = file.path(out_dir, "summary_input")
)
```

``` r

summary_fdr_preview <- data.table::fread(res_summary$file_fdr_out, nrows = 5)
for (column in c("best_hit", "best_hit_bonf", "best_hit_acat", "q")) {
  summary_fdr_preview[[column]] <- format(summary_fdr_preview[[column]], digits = 3, scientific = TRUE)
}
knitr::kable(
  summary_fdr_preview,
  caption = "Window-level results from summary-statistics input: first five rows."
)
```

| group | best_hit | n_var_in_cis | best_hit_bonf | best_hit_acat | snp          | q        |
|:------|:---------|-------------:|:--------------|:--------------|:-------------|:---------|
| G1    | 1.12e-02 |          563 | 1             | 9.92e-01      | 5:60134:T:A  | 9.98e-01 |
| G10   | 1.21e-02 |         1125 | 1             | 9.95e-01      | 5:632227:G:C | 9.98e-01 |
| G11   | 8.73e-04 |         1531 | 1             | 9.95e-01      | 5:945144:C:T | 9.98e-01 |
| G12   | 1.10e-02 |         1411 | 1             | 9.95e-01      | 5:853507:A:G | 9.98e-01 |
| G13   | 1.37e-02 |         1226 | 1             | 9.97e-01      | 5:975668:A:C | 9.98e-01 |

Window-level results from summary-statistics input: first five rows.
{.table}

The two modes return the same types of output files. These examples
demonstrate workflow execution; they do not evaluate statistical power,
calibration, or biological validity.

### Select chromosomes or skip FDR

The following recipe uses `chr = "chr5"` for a targeted analysis and
`do_fdr = FALSE` to return association P values without q-values. It is
not rerun here because the quick-start inputs already contain only
chromosome 5.

``` r

res_chr <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  qtl_file = qtl_file,
  chr = "chr5",
  do_fdr = FALSE,
  out_prefix = file.path(out_dir, "chr5_only")
)
```

For several chromosomes, supply a character vector to `chr`. `qtl_files`
accepts a combined summary-statistics file, one file per selected
chromosome, or a filename template containing `{chr}`. Use
[`?cacti_peak_window`](https://xuanyaoliulab.github.io/cacti/reference/cacti_peak_window.md),
[`?cacti_run_chr`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_chr.md),
and
[`?cacti_run_genome`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_genome.md)
for the full argument descriptions.
