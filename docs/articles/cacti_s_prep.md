# CACTI-S: segment-based data preparation & mapping

## Introduction

The **CACTI-S (Segment-based)** pipeline is designed to process raw
sequencing data (BAM files) into normalized phenotype matrices and QTL
mapping. It handles the multi-step pipeline:

- **Segmentation**: Tile the genome into fixed-size segments (e.g., 1kb,
  5kb).
- **Reads counting**: Count reads from raw BAM files into these
  segments.
- **Preprocessing**: Perform quality control and normalization.
- **Mapping**: Run cis-QTL mapping of the segments using `MatrixEQTL` to
  generate summary statistics required by the multivariate CACTI model.

This vignette explains:

1.  The overall workflow.
2.  The expected input files.
3.  The output files produced by the pipeline.
4.  Example usage with bundled toy data under `inst/extdata/`.

The main user-facing function in this vignette is:

- [`cacti_s_prep()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_prep.md)
  — a master wrapper that automates the entire workflow from BAM files
  to summary statistics.

Lower-level building blocks are:

| Function | Role | Description |
|:---|:---|:---|
| **`cacti_s_prep`** | **Wrapper** | Master function that runs the entire pipeline (Steps 1–4) sequentially. |
| `cacti_s_create_segments` | Step 1 | Generates fixed-size segments to divide the genome. |
| `cacti_s_count` | Step 2 | Quantifies reads in segments using [`Rsubread::featureCounts`](https://rdrr.io/pkg/Rsubread/man/featureCounts.html). |
| `cacti_s_preprocess` | Step 3 | Performs QC and normalization. |
| `cacti_s_map_cis` | Step 4 | Runs cis-QTL mapping using `MatrixEQTL`. |

## Part 1: Overview of the CACTI-S pipeline

The CACTI-S pipeline has three main steps:

1.  **Define and quantify genomic segments**

- Segmentation: Tiling the genome into fixed-size segments (e.g., 1kb,
  5kb) or using a custom BED file.

- Quantification: Counting reads in these segments for every sample (BAM
  file) using `Rsubread`.

- Output: A raw counts matrix (Rows = Segments, Columns = Samples).

2.  **Preprocess and normalize phenotypes**

- QC & Filtering & Normalization

- Output: A normalized phenotype matrix ready for QTL mapping.

3.  **Cis-QTLs mapping**

- Input: Normalized phenotype matrix, genotype data (VCF), and
  covariates.

- Process: Runs MatrixEQTL for every variant–segment pair within a
  defined cis window (e.g., 100kb).

- Output: Cis-QTL summary statistics (signed association statistics),
  which serve as the primary input for the CACTI multivariate test.

## Part 2: Quick start

For most users, the
[`cacti_s_prep()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_prep.md)
function is the most convenient way to run the pipeline. It automates
all steps and manages intermediate files.

### 1. Locate test data

We will use the small test files bundled with the `cacti` package. The
four BAM files contain distinct count profiles across twelve 5-kb
segments. These synthetic data demonstrate execution, not statistical
power or biological validity. In a real analysis, replace these with
your own BAM and VCF paths.

Identical warnings are displayed only once on this page. The generated
association statistics and P values are checked to be finite.

``` r

# Locate the 4 test BAM files
bam_files <- list.files(
  system.file("extdata", package = "cacti"),
  pattern = "Sample.*\\.bam$",
  full.names = TRUE
)

# Locate genotype and covariate files
file_vcf   <- system.file("extdata", "test_geno.vcf", package = "cacti")
file_cov   <- system.file("extdata", "test_cov.txt", package = "cacti")

print(basename(c(bam_files, file_vcf, file_cov)))
#> [1] "Sample1.bam"   "Sample2.bam"   "Sample3.bam"   "Sample4.bam"  
#> [5] "test_geno.vcf" "test_cov.txt"
```

### 2. Run the pipeline

We will run the pipeline with the following settings:

- **Genome**: hg19  
- **Segment Size**: 5kb  
- **Filter Mode**: `cpm_threshold` (keep segments with peak intensity
  CPM \> 1 in \> 20% of samples)

``` r

out_dir <- file.path(tempdir(), "cacti_s_example")
res <- cacti_s_prep(
  file_bams = bam_files,
  file_vcf = file_vcf,
  file_cov = file_cov,
  
  out_dir = out_dir,
  out_prefix = "test",

  # Pipeline Parameters
  genome = "hg19",
  segment_size = "5kb",
  filter_mode = "cpm_threshold",
  min_cpm = 1,
  min_prop = 0.2,

  # Use 1 thread for this small demo
  threads = 1
)
#> Created output directory: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example
#> =======================================================
#>  CACTI-S Preparation Pipeline
#> =======================================================
#> Output Directory: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example
#> Input BAMs: 4
#> 
#> [Step 1] Creating Genomic Segments...
#> Using built-in chromosome sizes for: hg19
#> Generating 5kb segments for 22 autosomes (chr1-chr22)...
#> Success! Created 576216 segments. Saved to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_segments.saf
#> 
#> [Step 2] Counting Reads...
#> Counting reads for 4 BAM files...
#> Settings: SAF format | Primary Only | Ignore Dups | Read2Pos = 5
#> 
#>         ==========     _____ _    _ ____  _____  ______          _____  
#>         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
#>           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
#>             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
#>               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
#>         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
#>        Rsubread 2.24.0
#> 
#> //========================== featureCounts setting ===========================\\
#> ||                                                                            ||
#> ||             Input files : 4 BAM files                                      ||
#> ||                                                                            ||
#> ||                           Sample1.bam                                      ||
#> ||                           Sample2.bam                                      ||
#> ||                           Sample3.bam                                      ||
#> ||                           Sample4.bam                                      ||
#> ||                                                                            ||
#> ||              Paired-end : no                                               ||
#> ||        Count read pairs : no                                               ||
#> ||              Annotation : test_segments.saf (SAF)                          ||
#> ||      Dir for temp files : .                                                ||
#> ||                 Threads : 1                                                ||
#> ||                   Level : meta-feature level                               ||
#> ||      Multimapping reads : counted                                          ||
#> ||     Multiple alignments : primary alignment only                           ||
#> || Multi-overlapping reads : not counted                                      ||
#> ||   Min overlapping bases : 1                                                ||
#> ||          Read reduction : to 5' end                                        ||
#> ||        Duplicated Reads : ignored                                          ||
#> ||                                                                            ||
#> \\============================================================================//
#> 
#> //================================= Running ==================================\\
#> ||                                                                            ||
#> || Load annotation file test_segments.saf ...                                 ||
#> ||    Features : 576216                                                       ||
#> ||    Meta-features : 576216                                                  ||
#> ||    Chromosomes/contigs : 22                                                ||
#> ||                                                                            ||
#> || Process BAM file Sample1.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 719                                                  ||
#> ||    Successfully assigned alignments : 719 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Process BAM file Sample2.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 674                                                  ||
#> ||    Successfully assigned alignments : 674 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Process BAM file Sample3.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 733                                                  ||
#> ||    Successfully assigned alignments : 733 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Process BAM file Sample4.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 590                                                  ||
#> ||    Successfully assigned alignments : 590 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Write the final count table.                                               ||
#> || Write the read assignment summary.                                         ||
#> ||                                                                            ||
#> \\============================================================================//
#> Success! Counts written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_raw_counts.txt
#> 
#> [Step 3] Preprocessing (QC, Filter, Normalize)...
#> === Starting Preprocessing ===
#>   Input: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_raw_counts.txt
#>   Filter Mode: cpm_threshold
#>   Cleaned sample names: Sample1, Sample2, Sample3...
#>   Initial Features: 576216
#>   Initial Samples:  4
#>   Features after removing all-zeros: 12
#>   Filtering by intensity: CPM > 1 in > 20% samples
#>   Features kept: 12
#>   Standardizing (Z-score per feature)...
#>   Rank Normalizing (INT per sample)...
#> Success! Processed data written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_pheno_norm.txt
#> Success! Processed data meta written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_pheno_norm_meta.txt
#> 
#> [Step 4] Mapping Cis-QTLs...
#> Preparing phenotype data...
#>   Phenotype samples: 4
#> Converting VCF to genotype matrix...
#> Preparing covariates...
#> Rows read: 51 done.
#> Rows read: 12 done.
#> Rows read: 1 done.
#> Running MatrixEQTL...
#> Matching data files and location files
#> 12 of 12 genes matched
#> 51 of 51 SNPs matched
#> Task finished in 0.001 seconds
#> Processing covariates
#> Task finished in 0 seconds
#> Processing gene expression data (imputation, residualization)
#> Task finished in 0.001 seconds
#> Creating output file(s)
#> Task finished in 0.001 seconds
#> Performing eQTL analysis
#> 100.00% done, 612 cis-eQTLs
#> Task finished in 0.004 seconds
#> 
#> Success! Stats written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_cis_qtl_stats.txt
#> Mapping complete: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_cis_qtl_stats.txt
#> 
#> =======================================================
#>  Pipeline Complete!
#> =======================================================
#> Summary of Generated Files:
#>   [1] Segments (SAF):      /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_segments.saf
#>   [2] Raw Counts:          /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_raw_counts.txt
#>   [3] Normalized Pheno:    /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_pheno_norm.txt
#>   [4] Cis-QTL Stats:       /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_cis_qtl_stats.txt
#> =======================================================
```

### 3. Explore results

The function returns a list of paths to the generated files. Outputs are
written to a temporary directory in this example; use a directory of
your choice to keep them. The table below shows the filenames.

``` r

knitr::kable(data.frame(Output = names(res), Filename = basename(unlist(res)), row.names = NULL))
```

| Output     | Filename               |
|:-----------|:-----------------------|
| saf        | test_segments.saf      |
| raw_counts | test_raw_counts.txt    |
| pheno      | test_pheno_norm.txt    |
| qtl_stats  | test_cis_qtl_stats.txt |

#### Segments

This file contains the segment positions that divide the whole genome.

``` r

saf <- read.table(res$saf, header = TRUE, nrows = 5)
knitr::kable(saf, caption = "Segment position (Top 5 Rows)")
```

| GeneID           | Chr  | Start |   End | Strand |
|:-----------------|:-----|------:|------:|:-------|
| chr1_1_5000      | chr1 |     1 |  5000 | .      |
| chr1_5001_10000  | chr1 |  5001 | 10000 | .      |
| chr1_10001_15000 | chr1 | 10001 | 15000 | .      |
| chr1_15001_20000 | chr1 | 15001 | 20000 | .      |
| chr1_20001_25000 | chr1 | 20001 | 25000 | .      |

Segment position (Top 5 Rows) {.table}

#### Raw counts

This file contains raw counts within the segments.

``` r

raw_counts <- data.table::fread(res$raw_counts, nrows = 5, data.table = FALSE)
names(raw_counts)[7:ncol(raw_counts)] <- basename(names(raw_counts)[7:ncol(raw_counts)])
knitr::kable(raw_counts, caption = "Raw counts within segments (Top 5 Rows)")
```

| GeneID | Chr | Start | End | Strand | Length | Sample1.bam | Sample2.bam | Sample3.bam | Sample4.bam |
|:---|:---|---:|---:|:---|---:|---:|---:|---:|---:|
| chr1_1_5000 | chr1 | 1 | 5000 | . | 5000 | 82 | 64 | 41 | 78 |
| chr1_5001_10000 | chr1 | 5001 | 10000 | . | 5000 | 89 | 51 | 86 | 36 |
| chr1_10001_15000 | chr1 | 10001 | 15000 | . | 5000 | 56 | 41 | 73 | 84 |
| chr1_15001_20000 | chr1 | 15001 | 20000 | . | 5000 | 10 | 83 | 63 | 30 |
| chr1_20001_25000 | chr1 | 20001 | 25000 | . | 5000 | 19 | 52 | 36 | 88 |

Raw counts within segments (Top 5 Rows) {.table}

#### Normalized phenotype matrix

This file contains filtered and normalized phenotypes, with segments in
rows and samples in columns. Segment positions are stored in the
separate metadata file.

``` r

pheno <- read.table(res$pheno, header = TRUE)
knitr::kable(head(pheno, 5), digits = 3, caption = "First five normalized segment phenotypes")
```

| GeneID           | Sample1 | Sample2 | Sample3 | Sample4 |
|:-----------------|--------:|--------:|--------:|--------:|
| chr1_1_5000      |  -0.105 |  -0.105 |  -1.732 |   0.812 |
| chr1_5001_10000  |   0.549 |  -0.549 |   1.150 |  -0.549 |
| chr1_10001_15000 |  -0.549 |  -1.150 |   0.319 |   1.732 |
| chr1_15001_20000 |  -1.150 |   1.150 |   0.549 |  -0.105 |
| chr1_20001_25000 |  -0.812 |   0.319 |  -0.812 |   1.150 |

First five normalized segment phenotypes {.table}

#### Cis-QTL summary statistics

This file contains the cis-QTL mapping associations of individual
segments, with signed Z-scores (`z`) and P values (`pval`).

``` r

if (!is.null(res$qtl_stats)) {
  qtl <- data.table::fread(res$qtl_stats, data.table = FALSE)
  qtl_preview <- head(qtl, 5)
  qtl_preview$pval <- formatC(qtl_preview$pval, format = "e", digits = 3)
  knitr::kable(qtl_preview, digits = 3, caption = "First five SNP-segment associations")
}
```

| phe_id           | var_id |      z | pval      |
|:-----------------|:-------|-------:|:----------|
| chr1_35001_40000 | snp_5  | -3.183 | 1.457e-03 |
| chr1_35001_40000 | snp_40 |  3.183 | 1.457e-03 |
| chr1_35001_40000 | snp_42 | -3.183 | 1.457e-03 |
| chr1_1_5000      | snp_13 |  2.700 | 6.931e-03 |
| chr1_1_5000      | snp_22 | -2.700 | 6.931e-03 |

First five SNP-segment associations {.table}

``` r

knitr::kable(data.frame(
  Output = c("Retained segments", "SNP-segment associations"),
  Number = c(nrow(pheno), nrow(qtl))
))
```

| Output                   | Number |
|:-------------------------|-------:|
| Retained segments        |     12 |
| SNP-segment associations |    612 |

## Part 3: Step-by-step workflow

If you need more control (e.g., inspecting intermediate files) or
customization, you can run the individual functions included in the
whole pipeline separately.

### Step 1: Create segments

First, we generate a SAF (Simplified Annotation Format) file that
defines the genomic segments.

``` r

file_saf <- file.path(out_dir, "test_segments.saf")

cacti_s_create_segments(
  file_saf_out = file_saf, # Output path for the SAF file.
  genome = "hg19", # Genome build to use: "hg19" (default) or "hg38".
  segment_size = "5kb" # Character like "5kb" or numeric bp (e.g. 5000).
)
#> Using built-in chromosome sizes for: hg19
#> Generating 5kb segments for 22 autosomes (chr1-chr22)...
#> Success! Created 576216 segments. Saved to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_segments.saf
```

``` r

head(read.table(file_saf, header = TRUE))
#>             GeneID  Chr Start   End Strand
#> 1      chr1_1_5000 chr1     1  5000      .
#> 2  chr1_5001_10000 chr1  5001 10000      .
#> 3 chr1_10001_15000 chr1 10001 15000      .
#> 4 chr1_15001_20000 chr1 15001 20000      .
#> 5 chr1_20001_25000 chr1 20001 25000      .
#> 6 chr1_25001_30000 chr1 25001 30000      .
```

### Step 2: Count Reads

Next, we quantify reads falling into these segments for all BAM files.

``` r

file_counts <- file.path(out_dir, "test_raw_counts.txt")

cacti_s_count(
  file_bams = bam_files, # Character vector of BAM file paths.
  file_saf = file_saf, # Path to the SAF annotation file (from \code{cacti_s_create_segments}).
  file_counts_out = file_counts, # Output path for the counts matrix.
  threads = 1 # Number of threads for featureCounts.
)
#> Counting reads for 4 BAM files...
#> Settings: SAF format | Primary Only | Ignore Dups | Read2Pos = 5
#> 
#>         ==========     _____ _    _ ____  _____  ______          _____  
#>         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
#>           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
#>             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
#>               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
#>         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
#>        Rsubread 2.24.0
#> 
#> //========================== featureCounts setting ===========================\\
#> ||                                                                            ||
#> ||             Input files : 4 BAM files                                      ||
#> ||                                                                            ||
#> ||                           Sample1.bam                                      ||
#> ||                           Sample2.bam                                      ||
#> ||                           Sample3.bam                                      ||
#> ||                           Sample4.bam                                      ||
#> ||                                                                            ||
#> ||              Paired-end : no                                               ||
#> ||        Count read pairs : no                                               ||
#> ||              Annotation : test_segments.saf (SAF)                          ||
#> ||      Dir for temp files : .                                                ||
#> ||                 Threads : 1                                                ||
#> ||                   Level : meta-feature level                               ||
#> ||      Multimapping reads : counted                                          ||
#> ||     Multiple alignments : primary alignment only                           ||
#> || Multi-overlapping reads : not counted                                      ||
#> ||   Min overlapping bases : 1                                                ||
#> ||          Read reduction : to 5' end                                        ||
#> ||        Duplicated Reads : ignored                                          ||
#> ||                                                                            ||
#> \\============================================================================//
#> 
#> //================================= Running ==================================\\
#> ||                                                                            ||
#> || Load annotation file test_segments.saf ...                                 ||
#> ||    Features : 576216                                                       ||
#> ||    Meta-features : 576216                                                  ||
#> ||    Chromosomes/contigs : 22                                                ||
#> ||                                                                            ||
#> || Process BAM file Sample1.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 719                                                  ||
#> ||    Successfully assigned alignments : 719 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Process BAM file Sample2.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 674                                                  ||
#> ||    Successfully assigned alignments : 674 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Process BAM file Sample3.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 733                                                  ||
#> ||    Successfully assigned alignments : 733 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Process BAM file Sample4.bam...                                            ||
#> ||    Single-end reads are included.                                          ||
#> ||    Total alignments : 590                                                  ||
#> ||    Successfully assigned alignments : 590 (100.0%)                         ||
#> ||    Running time : 0.00 minutes                                             ||
#> ||                                                                            ||
#> || Write the final count table.                                               ||
#> || Write the read assignment summary.                                         ||
#> ||                                                                            ||
#> \\============================================================================//
#> Success! Counts written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_raw_counts.txt
```

``` r

head(read.table(file_counts, header = TRUE))
#>             GeneID  Chr Start   End Strand Length Sample1.bam Sample2.bam
#> 1      chr1_1_5000 chr1     1  5000      .   5000          82          64
#> 2  chr1_5001_10000 chr1  5001 10000      .   5000          89          51
#> 3 chr1_10001_15000 chr1 10001 15000      .   5000          56          41
#> 4 chr1_15001_20000 chr1 15001 20000      .   5000          10          83
#> 5 chr1_20001_25000 chr1 20001 25000      .   5000          19          52
#> 6 chr1_25001_30000 chr1 25001 30000      .   5000          60          33
#>   Sample3.bam Sample4.bam
#> 1          41          78
#> 2          86          36
#> 3          73          84
#> 4          63          30
#> 5          36          88
#> 6          48          40
```

### Step 3: Preprocess

We then perform preprocessing steps on the count data, including -

- **QC:** Remove segments with 0 reads in all samples.  
- **CPM:** Transform counts to log2(CPM + 0.5).  
- **Filter:** Remove noise segments based on `filter_mode`. We use the
  default mode: `cpm_threshold` (keep segments with peak intensity CPM
  \> 1 in \> 20% of samples). We also provide other ways for filtering -
  - `match_overlap_count`: Determines $`N`$, the number of segments
    overlapping regions in `filter_bed`. Then keeps the top $`N`$
    segments by average intensity.
  - `prop`: Keeps the top `prop_top` proportion of segments by average
    intensity.
- **Standardize:** Z-score normalization per feature (Mean = 0, SD =
  1).  
- **RankNorm:** Inverse Normal Transform per sample.

``` r

file_pheno <- file.path(out_dir, "test_pheno_norm.txt")
file_pheno_meta <- file.path(out_dir, "test_pheno_norm_meta.txt")

cacti_s_preprocess(
  file_raw_counts = file_counts, # Path to featureCounts output (from \code{cacti_s_count}).
  file_out = file_pheno, # Output path for the processed phenotype file.
  file_out_meta = file_pheno_meta, # Output path for the meta file of processed phenotype.
  filter_mode = "cpm_threshold", # Strategy to select segments: "cpm_threshold" (default).
  min_cpm = 1, # (For "cpm_threshold") Minimum CPM threshold (default 1).
  min_prop = 0.2 # (For "cpm_threshold") Minimum proportion of samples (default 0.2).
)
#> === Starting Preprocessing ===
#>   Input: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_raw_counts.txt
#>   Filter Mode: cpm_threshold
#>   Cleaned sample names: Sample1, Sample2, Sample3...
#>   Initial Features: 576216
#>   Initial Samples:  4
#>   Features after removing all-zeros: 12
#>   Filtering by intensity: CPM > 1 in > 20% samples
#>   Features kept: 12
#>   Standardizing (Z-score per feature)...
#>   Rank Normalizing (INT per sample)...
#> Success! Processed data written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_pheno_norm.txt
#> Success! Processed data meta written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_pheno_norm_meta.txt
```

``` r

head(read.table(file_pheno, header = TRUE))
#>             GeneID    Sample1    Sample2    Sample3    Sample4
#> 1      chr1_1_5000 -0.1046335 -0.1046335 -1.7316644  0.8122178
#> 2  chr1_5001_10000  0.5485223 -0.5485223  1.1503494 -0.5485223
#> 3 chr1_10001_15000 -0.5485223 -1.1503494  0.3186394  1.7316644
#> 4 chr1_15001_20000 -1.1503494  1.1503494  0.5485223 -0.1046335
#> 5 chr1_20001_25000 -0.8122178  0.3186394 -0.8122178  1.1503494
#> 6 chr1_25001_30000  1.1503494 -1.7316644 -0.3186394  0.1046335
```

``` r

head(read.table(file_pheno_meta, header = TRUE))
#>             GeneID  Chr Start   End
#> 1      chr1_1_5000 chr1     1  5000
#> 2  chr1_5001_10000 chr1  5001 10000
#> 3 chr1_10001_15000 chr1 10001 15000
#> 4 chr1_15001_20000 chr1 15001 20000
#> 5 chr1_20001_25000 chr1 20001 25000
#> 6 chr1_25001_30000 chr1 25001 30000
```

### Step 4: Map Cis-QTLs

Finally, we run the cis-QTL mapping using the normalized phenotype
matrix using `MatrixEQTL`.

``` r

file_qtl <- file.path(out_dir, "test_cis_qtl_stats.txt")

cacti_s_map_cis(
  file_pheno = file_pheno, # Path to the processed phenotype file (from `cacti_s_preprocess`).
  file_pheno_meta = file_pheno_meta, # Path for the meta file of processed phenotype.
  file_cov = file_cov, # Path to covariate matrix.
  file_vcf = file_vcf, # Path to input VCF file.
  file_qtl_out = file_qtl, # Output path for summary statistics.
  cis_dist = 100000, # Cis-window distance (default 100000 bp = 100kb).
  p_threshold = 1.0 # P-value threshold for output (default 1.0, print all associations).
)
#> Preparing phenotype data...
#>   Phenotype samples: 4
#> Converting VCF to genotype matrix...
#> Preparing covariates...
#> Rows read: 51 done.
#> Rows read: 12 done.
#> Rows read: 1 done.
#> Running MatrixEQTL...
#> Matching data files and location files
#> 12 of 12 genes matched
#> 51 of 51 SNPs matched
#> Task finished in 0 seconds
#> Processing covariates
#> Task finished in 0.001 seconds
#> Processing gene expression data (imputation, residualization)
#> Task finished in 0 seconds
#> Creating output file(s)
#> Task finished in 0.002 seconds
#> Performing eQTL analysis
#> 100.00% done, 612 cis-eQTLs
#> Task finished in 0.005 seconds
#> 
#> Success! Stats written to: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_cis_qtl_stats.txt
```

``` r

head(read.table(file_qtl, header = TRUE))
#>             phe_id var_id         z        pval
#> 1 chr1_35001_40000  snp_5 -3.183188 0.001456629
#> 2 chr1_35001_40000 snp_40  3.183188 0.001456629
#> 3 chr1_35001_40000 snp_42 -3.183188 0.001456629
#> 4      chr1_1_5000 snp_13  2.700154 0.006930729
#> 5      chr1_1_5000 snp_22 -2.700154 0.006930729
#> 6      chr1_1_5000 snp_33  2.700154 0.006930729
```

## Part 4: Run multivariate association with CACTI

The files generated by the CACTI-S pipeline serve as the primary inputs
for the core CACTI analysis. These include the phenotype metadata
(`file_pheno_meta`), normalized phenotypes (`file_pheno`), covariates
(`file_cov`), and the QTL summary statistics (`file_qtl`). For more
details, please see the vignette
[here](https://xuanyaoliulab.github.io/cacti/articles/cacti_peak_window.md).

With these inputs ready, we can directly run
[`cacti_run_chr()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_chr.md)
(for a single chromosome) or
[`cacti_run_genome()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_genome.md)
(for the whole genome).

For the bundled example, 10-kb windows group adjacent 5-kb segments,
allowing windows with at least two retained segments to enter the
multivariate test.

``` r

file_pheno <- file.path(out_dir, "test_pheno_norm.txt")
file_pheno_meta <- file.path(out_dir, "test_pheno_norm_meta.txt")
file_cov <- system.file("extdata", "test_cov.txt", package = "cacti")
file_qtl <- file.path(out_dir, "test_cis_qtl_stats.txt")
out_prefix <- file.path(out_dir, "cactis_chr1")
chr <- "chr1"

res_single_chr <- cacti_run_chr(
  window_size     = "10kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno      = file_pheno,
  file_cov        = file_cov,
  chr             = chr,
  qtl_file        = file_qtl,
  out_prefix      = out_prefix,
  min_peaks       = 2
)
#> =======================================================
#>  CACTI Pipeline
#> =======================================================
#> Output Directory: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example
#> 
#> [Step 1/3] Group peaks into non-overlapping windows …
#> 
#> [Step 2/3] Residualize phenotype by covariates …
#> 
#> [Step 3/3] Run pvalue for chr1 …
#> Warning in ACAT::ACAT(rbind(p.PCMinP, p.PCFisher, p.PCLC, p.WI, p.Wald, : There
#> are p-values that are exactly 1!
#> 
#> =======================================================
#> 
#>  CACTI pipeline for chr1 completed!
#> =======================================================
#> =======================================================
#> ----------- Run summary -----------
#> Chromosome:chr1
#> Window size:10kb
#> 
#>  ----------- Input files -----------
#> Peaks BED:/var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_pheno_norm_meta.txt
#> Phenotype:/var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_pheno_norm.txt
#> Covariate:/private/var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T/RtmpO6YTYF/temp_libpath11ca5d430a3c/cacti/extdata/test_cov.txt
#> QTL summary stats:/var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/test_cis_qtl_stats.txt
#> 
#>  ----------- Output files -----------
#>   [1] Grouped peak: /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/cactis_chr1_peak_group_window10kb.txt
#>   [2] Grouped peak (peak as row):/var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/cactis_chr1_peak_group_window10kb_peak_as_row.txt
#>   [3] Residual phenotype:/var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/cactis_chr1_pheno_cov_residual.txt
#>   [4] Pval results (chr1): /var/folders/p0/bl6jnkv13_s1l9r9qy4dsndc0000gn/T//RtmpLPnDhi/cacti_s_example/cactis_chr1_pval_window10kb_chr1.txt.gz
#> =======================================================
```

The function returns a list containing the paths to the generated files:

``` r

knitr::kable(data.frame(Output = names(res_single_chr), Filename = basename(unlist(res_single_chr)), row.names = NULL))
```

| Output | Filename |
|:---|:---|
| file_peak_group | cactis_chr1_peak_group_window10kb.txt |
| file_peak_group_peaklevel | cactis_chr1_peak_group_window10kb_peak_as_row.txt |
| file_pheno_cov_residual | cactis_chr1_pheno_cov_residual.txt |
| file_p_peak_group | cactis_chr1_pval_window10kb_chr1.txt.gz |

Below is a preview of the final association statistics:

``` r

window_results <- data.table::fread(res_single_chr$file_p_peak_group, data.table = FALSE)
window_preview <- head(window_results, 5)
window_preview$pval <- formatC(window_preview$pval, format = "e", digits = 3)
knitr::kable(window_preview, caption = "First five SNP-window associations")
```

| group | snp    | pval      |
|:------|:-------|:----------|
| G1    | snp_13 | 3.056e-02 |
| G1    | snp_22 | 3.056e-02 |
| G1    | snp_33 | 3.056e-02 |
| G1    | snp_44 | 3.056e-02 |
| G1    | snp_11 | 3.056e-02 |

First five SNP-window associations {.table}

``` r

knitr::kable(data.frame(
  Output = c("Tested windows", "SNP-window associations"),
  Number = c(length(unique(window_results$group)), nrow(window_results))
))
```

| Output                  | Number |
|:------------------------|-------:|
| Tested windows          |      6 |
| SNP-window associations |    306 |
