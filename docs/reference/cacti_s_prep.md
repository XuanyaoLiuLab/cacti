# Run CACTI-S Prep Pipeline

A master wrapper function that executes the full CACTI-S preparation
pipeline. It automates the transition from raw BAM files to a normalized
phenotype matrix and optionally performs nominal cis-QTL mapping.

## Usage

``` r
cacti_s_prep(
  file_bams,
  file_vcf = NULL,
  file_cov = NULL,
  out_dir,
  out_prefix = "cacti_s",
  file_saf_custom = NULL,
  genome = "hg19",
  segment_size = "5kb",
  filter_mode = "cpm_threshold",
  filter_bed = NULL,
  min_cpm = 1,
  min_prop = 0.2,
  prop_top = 0.2,
  cis_dist = 1e+05,
  threads = 1,
  isPairedEnd = FALSE
)
```

## Arguments

- file_bams:

  Character vector of file paths to the input BAM files. These should be
  sorted and indexed.

- file_vcf:

  (Optional) Path to the input VCF file containing genotype data.
  Required if you want to run the mapping step.

- file_cov:

  (Optional) Path to the covariate matrix file. Format: Rows are
  covariates (ID, S1, S2...), Columns are samples. Required if
  `file_vcf` is provided.

- out_dir:

  Directory path where all output files will be saved. Created
  automatically if it does not exist.

- out_prefix:

  String prefix for all generated filenames (default "cacti_s"). Files
  created: `_segments.saf`, `_raw_counts.txt`, `_pheno_norm.txt`,
  `_cis_qtl_stats.txt`, etc.

- file_saf_custom:

  (Optional) Path to an existing SAF file. If provided, segment
  generation is skipped.

- genome:

  Genome build to use for segment generation. Options: `"hg19"`
  (default) or `"hg38"`.

- segment_size:

  String defining the window size for tiling the genome (e.g., `"1kb"`,
  `"5kb"`, `"10kb"`). Default is `"5kb"`.

- filter_mode:

  Strategy for selecting informative segments during preprocessing:

  - `"cpm_threshold"` (Default): Keeps segments with Counts Per Million
    (CPM) \> `min_cpm` in at least `min_prop` proportion of samples.

  - `"match_overlap_count"`: Calculates \$N\$, the number of segments
    that physically overlap with regions in `filter_bed`. It then
    selects the top \$N\$ segments based on average CPM intensity.

  - `"prop"`: Simply keeps the top `prop_top` proportion of segments
    ranked by average CPM.

- filter_bed:

  (Required only for `match_overlap_count`) Path to a BED file (e.g.,
  ATAC-seq peaks, H3K27ac peaks, or known QTL regions) used to determine
  the number of segments to keep.

- min_cpm:

  (For `cpm_threshold`) Minimum CPM intensity threshold. Default 1.

- min_prop:

  (For `cpm_threshold`) Minimum proportion of samples (0-1) that must
  exceed `min_cpm`. Default 0.2.

- prop_top:

  (For `prop`) Proportion (0-1) of top-expressed segments to keep.
  Default 0.2.

- cis_dist:

  (For Mapping) The maximum distance (bp) between a gene (segment) and a
  SNP to be considered "cis". Default 100,000 bp.

- threads:

  Integer, number of threads to use for the read counting step. Default
  1.

- isPairedEnd:

  Logical, indicating if the BAM files contain paired-end reads. Default
  `FALSE`.

## Value

Invisibly returns a list containing paths to the generated files:

- `segments.saf`: Path to the segment definition file.

- `raw_counts.txt`: Path to the featureCounts output.

- `pheno_norm.txt`: Path to the final normalized phenotype matrix.

- `cis_qtl_stats.txt`: (If mapped) Path to the cis-QTL summary
  statistics.

## Details

**Pipeline Stages:**

1.  **Create Segments:** Generates fixed-size genomic windows (SAF) for
    the specified genome build.

2.  **Count Reads:** Quantifies reads in these segments using
    [`Rsubread::featureCounts`](https://rdrr.io/pkg/Rsubread/man/featureCounts.html).

3.  **Preprocess:** Performs QC (removing zero-count segments),
    Normalization (CPM), Filtering (different modes), Standardization
    (Z-score), and Rank Normalization (INT).

4.  **Map Cis-QTLs:** (Optional) Runs linear regression using
    `MatrixEQTL` if VCF and Covariates are provided.

## Examples

``` r
if (FALSE) { # \dontrun{
# Locate test data bundled with the package
file_bam1 <- system.file("extdata", "Sample1.bam", package = "cacti")
file_bam2 <- system.file("extdata", "Sample2.bam", package = "cacti")
file_bam3 <- system.file("extdata", "Sample3.bam", package = "cacti")
file_bam4 <- system.file("extdata", "Sample4.bam", package = "cacti")
file_vcf  <- system.file("extdata", "test_geno.vcf", package = "cacti")
file_cov  <- system.file("extdata", "test_cov.txt", package = "cacti")

# Create a temporary output directory
out_dir <- "inst/extdata/test_results"

if (nchar(file_bam1) > 0) {
  # --- Run Full Pipeline (BAM -> QTL Stats) ---
  res <- cacti_s_prep(
    file_bams = c(file_bam1, file_bam2, file_bam3, file_bam4),
    file_vcf = file_vcf,
    file_cov = file_cov,

    out_dir = out_dir,
    out_prefix = "test",

    # Segment Parameters
    file_saf_custom = NULL,
    genome = "hg19",
    segment_size = "5kb",

    # Filtering Parameters
    filter_mode = "cpm_threshold",
    min_cpm = 1,
    min_prop = 0.2,

    # Execution Parameters
    threads = 1
  )
}
} # }
```
