
# CACTI: Leveraging Correlated Regulatory Elements for Powerful Chromatin QTL Detection

## Overview

<img src="man/figures/cacti_overview.png" width="2055" />

CACTI
implements a powerful method for chromatin QTL mapping that leverages the correlation structure of nearby regulatory elements.


The package offers two main modules:
    
1.  **CACTI (Peak-based pipeline):**
* Suitable for chromatin features with narrow peak signatures (H3k27ac, H3k4me1, ATAC-seq, etc)
* Peak calling should be performed before applying CACTI.
* Takes standard input in standard QTL calling tools: genotype, phenotype (called peaks) and covariates
* A summary statistics-based option is available too

2.  **CACTI-S (Segment-based pipeline):**
* Suitable for chromatin features with broad peak signatures (H3K27me3, H3K36me3, etc)
* Skips peak calling, an end-to-end preprocessing workflow to go from raw BAM files to normalized phenotype matrices to QTL mapping.
* Performs segmentation, read counting, QC/filtering, and normalization.
  

------------------------------------------------------------------------

## Installation

### Prerequisites

Before installing CACTI, ensure you have:

- **R >= 4.0.0** - Check your version with `R.version.string`
- **Pandoc** - Required for building vignettes. See [pandoc github](https://github.com/cderv/pandoc) for detailed instruction. 
- **ACAT** - Required for association testing. See [ACAT github](https://github.com/yaowuliu/ACAT) for detailed instruction.
- **qvalue** - Required for FDR correction. See [qvalue github](https://github.com/StoreyLab/qvalue) for detailed instruction.
- **Rsubread** - Required for read counting in CACTI-S. See [Rsubread website](https://www.bioconductor.org/packages//2.10/bioc/html/Rsubread.html) for detailed instruction.



``` r
install.packages("devtools")   # if not installed yet
devtools::install_github("XuanyaoLiuLab/cacti", build_vignettes = TRUE)
```

After installation:

``` r
library(cacti)
```

------------------------------------------------------------------------

## Getting Started in 5 Minutes

### 1) Run CACTI peak-window (genome-wide default)

```r
library(cacti)

file_pheno_meta <- system.file("extdata", "test_cacti_peak_chr5_pheno_meta.bed", package = "cacti")
file_pheno <- system.file("extdata", "test_cacti_peak_chr5_pheno.txt", package = "cacti")
file_cov <- system.file("extdata", "test_cacti_peak_chr5_covariates.txt", package = "cacti")
file_vcf <- system.file("extdata", "test_cacti_peak_chr5_geno.vcf", package = "cacti")
qtl_file <- system.file("extdata", "test_cacti_peak_chr5_matrixqtl_sumstats.txt.gz", package = "cacti")

res <- cacti_peak_window(
  window_size = "50kb", # window size to group peaks
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  file_vcf = file_vcf,
  chr = "All",      # default: run all chromosomes in file_pheno_meta
  do_fdr = TRUE,    # default: add FDR correction
  out_prefix = tempfile("cacti_quickstart_")
)

res$file_fdr_out

# Optional summary-stats mode
res_qtl <- cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno = file_pheno,
  file_cov = file_cov,
  qtl_file = qtl_file,
  chr = "chr5",
  do_fdr = TRUE,
  out_prefix = tempfile("cacti_quickstart_qtl_")
)

res_qtl$file_fdr_out
```

### 2) Run CACTI-S prep (raw BAM to cis-QTL summary stats)

```r
library(cacti)

file_bams <- c(
  system.file("extdata", "Sample1.bam", package = "cacti"),
  system.file("extdata", "Sample2.bam", package = "cacti"),
  system.file("extdata", "Sample3.bam", package = "cacti"),
  system.file("extdata", "Sample4.bam", package = "cacti")
)
file_vcf <- system.file("extdata", "test_geno.vcf", package = "cacti")
file_cov <- system.file("extdata", "test_cov.txt", package = "cacti")

res_s <- cacti_s_prep(
  file_bams = file_bams,
  file_vcf = file_vcf,
  file_cov = file_cov,
  out_dir = tempdir(),
  out_prefix = "cacti_s_quickstart"
)

res_s
```

Notes:
- Input conventions are standard QTL inputs: genotype, phenotype, and covariates with matched sample IDs.
- `cacti_peak_window()` also supports summary-statistics mode via `qtl_file`/`qtl_files`.

------------------------------------------------------------------------

## Documentation

### Vignettes

See the [full documentation and vignettes](https://xuanyaoliulab.github.io/cacti/) for -

- CACTI Peak-Window Pipeline

- CACTI-S Pipeline

Vignettes can also be viewed within the installed package -

```r
vignette("cacti_peak_window", package = "cacti")
vignette("cacti_s_prep", package = "cacti")
```


### Main functions documentation

```r
?cacti_run_chr
?cacti_run_genome
?cacti_add_fdr
?cacti_s_prep
```

------------------------------------------------------------------------


## Citation

If you use the CACTI method, please cite:

> Wang, L., & Liu, X. (2025). Improved chromatin QTL mapping with CACTI.
> bioRxiv, 2025-06.
