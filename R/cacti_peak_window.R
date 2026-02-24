#' Run CACTI peak-window pipeline end-to-end (genotype-first by default)
#'
#' High-level wrapper for the CACTI peak-window workflow.
#' By default, this function runs MatrixEQTL first from genotype + phenotype +
#' covariates, then runs standard CACTI peak-window testing.
#'
#' Summary-statistics mode is optional: if `qtl_files` (or `qtl_file`) is
#' provided, MatrixEQTL is skipped and those summary stats are used directly.
#'
#' @param window_size Character like `"50kb"` (kb units) or numeric (bp).
#' @param file_pheno_meta Path to phenotype meta table with columns:
#'   `phe_id`, `phe_chr`, `phe_from`, `phe_to`.
#' @param file_pheno Path to phenotype matrix. Rows = features, first column = ID.
#' @param file_cov Path to covariate matrix. Rows = covariates, first column = ID.
#' @param out_prefix Output prefix.
#' @param chrs Optional chromosome vector. If `NULL`, inferred from `file_pheno_meta`.
#' @param qtl_files Optional summary-stat input (vector, single path, or `{chr}` template).
#'   If `NULL`, MatrixEQTL is run first.
#' @param qtl_file Optional single-chromosome alias of `qtl_files`.
#' @param file_vcf Optional VCF genotype input for MatrixEQTL-first mode.
#' @param file_geno Optional genotype matrix input (if no VCF).
#' @param file_snp_pos Optional SNP-position input (if no VCF).
#' @param cis_dist Cis distance for MatrixEQTL (default 100kb).
#' @param p_threshold P-value threshold for MatrixEQTL output (default 1.0).
#' @param dir_pco Path to PCO helper files.
#' @param min_peaks Minimum number of peaks required for a window to be included in testing.
#' Included windows with 1 peak use univariate p-values; included windows with >=2 peaks use PCO.
#' @param file_fdr_out Optional output path for FDR-added result file.
#'
#' @return Invisibly returns the output list from [cacti_run_genome()].
#' @export
cacti_peak_window <- function(
    window_size = "50kb",
    file_pheno_meta,
    file_pheno,
    file_cov,
    out_prefix,
    chrs = NULL,
    qtl_files = NULL,
    qtl_file = NULL,
    file_vcf = NULL,
    file_geno = NULL,
    file_snp_pos = NULL,
    cis_dist = 100000,
    p_threshold = 1.0,
    dir_pco = system.file("pco", package = "cacti"),
    min_peaks = 1,
    file_fdr_out = NULL
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required.")

  if (!is.null(qtl_file) && !is.null(qtl_files)) {
    stop("Provide only one of `qtl_file` or `qtl_files`.")
  }
  if (!is.null(qtl_file)) qtl_files <- qtl_file

  if (is.null(chrs)) {
    meta_chr <- data.table::fread(file_pheno_meta, header = TRUE, select = "phe_chr")
    chrs <- unique(meta_chr$phe_chr)
  }

  if (length(chrs) == 0) stop("No chromosomes found in `file_pheno_meta`.")

  if (is.null(qtl_files) && is.null(file_vcf) && is.null(file_geno)) {
    stop(
      "Genotype-first mode requires `file_vcf` or `file_geno` when `qtl_files` is NULL. ",
      "Alternatively provide `qtl_files`/`qtl_file`."
    )
  }

  cacti_run_genome(
    window_size = window_size,
    file_pheno_meta = file_pheno_meta,
    file_pheno = file_pheno,
    file_cov = file_cov,
    chrs = chrs,
    qtl_files = qtl_files,
    file_vcf = file_vcf,
    file_geno = file_geno,
    file_snp_pos = file_snp_pos,
    cis_dist = cis_dist,
    p_threshold = p_threshold,
    out_prefix = out_prefix,
    dir_pco = dir_pco,
    min_peaks = min_peaks,
    file_fdr_out = file_fdr_out
  )
}
