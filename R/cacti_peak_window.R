#' Run CACTI peak-window pipeline end-to-end
#'
#' High-level wrapper for the CACTI peak-window workflow.
#' This function supports two input modes:
#' 1) genotype + phenotype + covariates (runs MatrixEQTL first),
#' 2) precomputed summary statistics (skips MatrixEQTL).
#'
#' By default (`chr = "All"`), it runs the genome-wide CACTI peak-window
#' workflow via [cacti_run_genome()] across all chromosomes present in
#' `file_pheno_meta` and performs FDR correction.
#'
#' @param window_size Size of non-overlapping genomic windows used to group
#'   peaks before CACTI testing. Accepts values like `"50kb"` (kilobases) or
#'   numeric base-pair values (e.g., `50000`).
#' @param file_pheno_meta Path to phenotype meta table with columns:
#'   `phe_id`, `phe_chr`, `phe_from`, `phe_to`.
#' @param file_pheno Path to phenotype matrix. Rows = features, first column = ID.
#' @param file_cov Path to covariate matrix. Rows = covariates, first column = ID.
#' @param out_prefix Output prefix.
#' @param chr Chromosome selection. Use `"All"` (default) to iterate over all
#'   chromosomes in `file_pheno_meta`, or provide one chromosome label
#'   (e.g., `"chr5"`), or a character vector of labels.
#' @param qtl_files Optional summary-stat input (vector, single path, or `{chr}` template).
#'   If `NULL`, MatrixEQTL is run first from genotype input.
#' @param qtl_file Optional single-chromosome alias of `qtl_files`.
#' @param file_vcf Optional VCF genotype input for genotype-input mode.
#' @param file_geno Optional genotype matrix input for genotype-input mode (if no VCF).
#' @param file_snp_pos Optional SNP-position input for genotype-input mode (if no VCF).
#' @param cis_dist Cis distance for MatrixEQTL (default 100kb).
#' @param dir_pco Path to PCO helper files.
#' @param min_peaks Minimum number of peaks required for a window to be included in testing.
#' Included windows with 1 peak use univariate p-values; included windows with >=2 peaks use PCO.
#' @param do_fdr Logical; if `TRUE` (default), run FDR correction across all
#'   processed chromosomes. If `FALSE`, skip FDR correction.
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
    chr = "All",
    qtl_files = NULL,
    qtl_file = NULL,
    file_vcf = NULL,
    file_geno = NULL,
    file_snp_pos = NULL,
    cis_dist = 100000,
    dir_pco = system.file("pco", package = "cacti"),
    min_peaks = 1,
    do_fdr = TRUE,
    file_fdr_out = NULL
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required.")

  if (!is.null(qtl_file) && !is.null(qtl_files)) {
    stop("Provide only one of `qtl_file` or `qtl_files`.")
  }
  if (!is.null(qtl_file)) qtl_files <- qtl_file

  meta_chr <- data.table::fread(file_pheno_meta, header = TRUE, select = "phe_chr")
  all_chrs <- unique(meta_chr$phe_chr)
  if (length(all_chrs) == 0) stop("No chromosomes found in `file_pheno_meta`.")

  if (length(chr) == 1L && is.character(chr) && tolower(chr) == "all") {
    chrs <- all_chrs
  } else {
    chrs <- as.character(chr)
  }
  if (length(chrs) == 0) stop("No chromosomes selected.")

  missing_chrs <- setdiff(chrs, all_chrs)
  if (length(missing_chrs) > 0) {
    stop("These chromosomes are not present in `file_pheno_meta`: ",
         paste(missing_chrs, collapse = ", "))
  }

  if (is.null(qtl_files) && is.null(file_vcf) && is.null(file_geno)) {
    stop(
      "Genotype-input mode requires `file_vcf` or `file_geno` when `qtl_files` is NULL. ",
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
    p_threshold = 1.0,
    out_prefix = out_prefix,
    dir_pco = dir_pco,
    min_peaks = min_peaks,
    file_fdr_out = file_fdr_out,
    do_fdr = do_fdr
  )
}
