# Create fixed-size genomic segments (SAF format)

Splits autosomes (chr1-chr22) into fixed-size segments and writes a SAF
file compatible with
[`Rsubread::featureCounts`](https://rdrr.io/pkg/Rsubread/man/featureCounts.html).

## Usage

``` r
cacti_s_create_segments(
  file_saf_out,
  genome = c("hg19", "hg38"),
  segment_size = "5kb"
)
```

## Arguments

- file_saf_out:

  Output path for the SAF file.

- genome:

  Genome build to use: "hg19" (default) or "hg38".

- segment_size:

  Character like "5kb" or numeric bp (e.g. 5000). Default is "5kb".

## Value

Invisibly returns the segment data frame. Writes `file_saf_out` to disk.

## Details

This function uses built-in chromosome sizes for hg19 or hg38. It
strictly filters the genome to include **only** chromosomes named "chr1"
through "chr22". Sex chromosomes (chrX, chrY) and mitochondrial DNA
(chrM) are excluded.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create 5kb segments for hg19 (default)
cacti_s_create_segments("inst/extdata/test_results/test_segments.saf")

# Create 1kb segments for hg38
cacti_s_create_segments("inst/extdata/test_results/test_segments.saf", genome = "hg38", segment_size = "1kb")
} # }
```
