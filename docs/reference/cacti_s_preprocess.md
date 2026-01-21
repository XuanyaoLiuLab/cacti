# Preprocess CACTI-S Segments

Performs QC, filtering, and normalization

## Usage

``` r
cacti_s_preprocess(
  file_raw_counts,
  file_out,
  file_out_meta,
  filter_mode = c("cpm_threshold", "match_overlap_count", "prop"),
  filter_bed = NULL,
  prop_top = 0.2,
  min_cpm = 1,
  min_prop = 0.2
)
```

## Arguments

- file_raw_counts:

  Path to featureCounts output (from `cacti_s_count`).

- file_out:

  Output path for the processed phenotype file.

- file_out_meta:

  Output path for the meta file of processed phenotype.

- filter_mode:

  Strategy to select segments: "cpm_threshold" (default),
  "match_overlap_count", or "prop".

- filter_bed:

  (Required for "match_overlap_count") Path to a BED file (e.g. peaks or
  QTL regions).

- prop_top:

  (For "prop") Proportion of segments to keep (0.0 to 1.0). Default 0.2.

- min_cpm:

  (For "cpm_threshold") Minimum CPM threshold (default 1).

- min_prop:

  (For "cpm_threshold") Minimum proportion of samples (default 0.2).

## Value

Invisibly returns the processed data frame. Writes `file_out` to disk.

## Details

**Processing Steps:**

1.  **QC:** Remove segments with 0 reads in all samples.

2.  **CPM:** Transform counts to log2(CPM + 0.5).

3.  **Filter:** Subset segments based on `filter_mode`.

4.  **Standardize:** Z-score normalization per feature (Mean=0, SD=1).

5.  **RankNorm:** Inverse Normal Transform per sample.

**Filtering Logic:**

- **cpm_threshold** (Default): Keeps segments with peak intensity \>
  `min_cpm` in \> `min_prop` of samples. Standard method to remove
  noise.

- **match_overlap_count**: Determines \$N\$, the number of segments
  overlapping regions in `filter_bed`. Then keeps the top \$N\$ segments
  by average intensity.

- **prop**: Keeps the top `prop_top` proportion of segments by average
  intensity.

## Examples

``` r
if (FALSE) { # \dontrun{
# Locate the bundled test data
filecounts <- system.file("extdata/test_results", "test_raw_counts.txt", package = "cacti")
filepeaks  <- system.file("extdata", "test_peaks.bed", package = "cacti")

# Define output file
fileout <- "inst/extdata/test_results/test_pheno_norm.txt"
fileout <- "inst/extdata/test_results/test_pheno_norm_meta.txt"

# Check if files exist (only true if package is installed)
if (nchar(filecounts) > 0) {

  # --- Example 1: Default Filtering (CPM Threshold) ---
  cacti_s_preprocess(
    file_raw_counts = filecounts,
    file_out = fileout,
    file_out_meta = fileout_meta,
    filter_mode = "cpm_threshold",
    min_cpm = 10,
    min_prop = 0.5
  )
  head(read.table(fileout, header = TRUE))

  # --- Example 2: Match Overlap Count ---
  # Logic: If 2 segments overlap peaks, keep the top 2 segments by intensity.
  cacti_s_preprocess(
    file_raw_counts = filecounts,
    file_out = fileout,
    file_out_meta = fileout_meta,
    filter_mode = "match_overlap_count",
    filter_bed = filepeaks
  )
  head(read.table(fileout, header = TRUE))
}
} # }
```
