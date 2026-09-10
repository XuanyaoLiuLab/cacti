# Package index

## CACTI peak-window cQTL mapping pipeline

- [`cacti_peak_window()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_peak_window.md)
  : Run CACTI peak-window pipeline end-to-end
- [`cacti_run_chr()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_chr.md)
  : Run CACTI peak-window pipeline for a single chromosome
- [`cacti_run_genome()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_run_genome.md)
  : Run CACTI peak-window pipeline genome-wide

## CACTI-S (Segment-based data preparation & mapping)

- [`cacti_s_prep()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_prep.md)
  : Run CACTI-S Prep Pipeline

## Building blocks

Low-level functions for segmentation, counting, preprocessing, grouping,
and testing.

- [`cacti_group_peak_window()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_group_peak_window.md)
  : Group adjacent peaks by non-overlapping windows
- [`cacti_pheno_cov_residual()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_pheno_cov_residual.md)
  : Residualize peak intensity with covariates
- [`cacti_cal_p()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_cal_p.md)
  : Run CACTI association test for a single chromosome
- [`cacti_add_fdr()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_add_fdr.md)
  : Add FDR (q-values) to per-window top hits across chromosomes
- [`cacti_s_create_segments()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_create_segments.md)
  : Create fixed-size genomic segments (SAF format)
- [`cacti_s_count()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_count.md)
  : Count reads in genomic segments (CACTI-S)
- [`cacti_s_preprocess()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_preprocess.md)
  : Preprocess CACTI-S Segments
- [`cacti_s_map_cis()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_s_map_cis.md)
  : Run Nominal Cis-QTL Mapping (CACTI-S)
- [`cacti_matrixqtl_cis()`](https://xuanyaoliulab.github.io/cacti/reference/cacti_matrixqtl_cis.md)
  : Run MatrixEQTL and write CACTI-compatible cis-QTL summary stats
