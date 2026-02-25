# CACTI: Chromatin QTL Mapping with Peak Windows

The `cacti` package implements the CACTI peak-window method for
chromatin QTL mapping. It is designed to:

## Details

- group nearby peaks into non-overlapping genomic windows,

- residualize peak intensities with respect to covariates,

- run multi-peak association tests using a principal-component omnibus
  (PCO) test, and

- aggregate per-chromosome results and compute window-level FDR.

The main user-facing functions are:

- [`cacti_run_chr`](https://liliw-w.github.io/cacti/reference/cacti_run_chr.md):
  run the peak-window pipeline for a single chromosome.

- [`cacti_run_genome`](https://liliw-w.github.io/cacti/reference/cacti_run_genome.md):
  run the peak-window pipeline across multiple chromosomes and add FDR.

Lower-level building blocks:

- [`cacti_group_peak_window`](https://liliw-w.github.io/cacti/reference/cacti_group_peak_window.md):
  group peaks into non-overlapping windows.

- [`cacti_pheno_cov_residual`](https://liliw-w.github.io/cacti/reference/cacti_pheno_cov_residual.md):
  residualize phenotypes by covariates.

- [`cacti_cal_p`](https://liliw-w.github.io/cacti/reference/cacti_cal_p.md):
  run CACTI association tests on peak windows.

- [`cacti_add_fdr`](https://liliw-w.github.io/cacti/reference/cacti_add_fdr.md):
  compute window-level FDR across chromosomes.

A small toy dataset from Aracena et al. is bundled under `inst/extdata/`
and can be accessed using
[`system.file`](https://rdrr.io/r/base/system.file.html) for testing and
examples.

For an overview of the peak-window pipeline and file formats, see the
package vignette:

`browseVignettes("cacti")`

## Author

**Maintainer**: Lili Wang <liliw1@uchicago.edu>

**Contributors**: Lili Wang, Zining Qi, Xuanyao Liu
