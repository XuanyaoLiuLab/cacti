#' Run CACTI peak-window pipeline for a single chromosome
#'
#' Steps:
#' 1) Group peaks into non-overlapping windows.
#' 2) Residualize phenotype by covariates.
#' 3) Run association test on peak groups for a single chromosome.
#'
#' @param window_size Character like `"50kb"` (kb units) or numeric (bp).
#' @param file_pheno_meta Path to input phenotype meta table (with header). Columns are -
#' - `phe_chr`: chromosome, with "chr" prefix (e.g. "chr1")
#' - `phe_from`: peak start (integer, < `phe_to`)
#' - `phe_to`: peak end (integer)
#' - `phe_id`: peak identifier (unique per row)
#' @param file_pheno Path to phenotype expression matrix.
#' Rows = features; first column = feature ID; remaining columns = samples.
#' @param file_cov Path to covariate matrix.
#' Rows = covariates; first column = covariate ID; remaining columns = samples.
#' @param chr Chromosome label to process (e.g., "chr1", "chr21").
#' @param qtl_file Optional path to cis QTL summary stats file for this chromosome.
#' Must contain columns: phe_id, var_id, z. If `NULL`, MatrixEQTL is run first
#' from genotype + phenotype + covariates to generate this file.
#' @param file_vcf Optional path to input VCF file for MatrixEQTL preprocessing.
#' @param file_geno Optional path to genotype matrix if no VCF is provided.
#' @param file_snp_pos Optional path to SNP-position file if no VCF is provided.
#' @param cis_dist Cis-window distance for MatrixEQTL (default 100000 bp).
#' @param p_threshold P-value threshold for MatrixEQTL output (default 1.0).
#' @param out_prefix Output prefix used to construct all output filenames.
#' @param dir_pco Directory containing association test helpers: ModifiedPCOMerged_acat.R, liu.R, liumod.R, davies.R, qfc.so.
#' @param min_peaks Minimum number of peaks required for a window to be included in testing.
#' Included windows with 1 peak use univariate p-values; included windows with >=2 peaks use PCO.
#'
#' @return Invisibly returns a named list of output paths with elements:
#'   \describe{
#'     \item{file_peak_group}{Path to the window-level group.}
#'     \item{file_peak_group_peaklevel}{Path to the peak-level group.}
#'     \item{file_pheno_residual}{Path to the residualized phenotype matrix.}
#'     \item{file_p_peak_group}{Path to the per-window p-value file for this chromosome.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example (chr5 only)
#'
#' file_pheno_meta <- system.file(
#'   "extdata", "test_pheno_meta.bed",
#'   package = "cacti"
#' )
#'
#' file_pheno <- system.file(
#'   "extdata", "test_pheno.txt",
#'   package = "cacti"
#' )
#'
#' file_cov <- system.file(
#'   "extdata", "test_covariates.txt",
#'   package = "cacti"
#' )
#'
#' out_prefix <- tempfile("cacti_chr5_")
#'
#' res <- cacti_run_chr(
#'   window_size = "50kb",
#'   file_pheno_meta = file_pheno_meta,
#'   file_pheno = file_pheno,
#'   file_cov = file_cov,
#'   chr = "chr5",
#'   qtl_file = system.file("extdata", "test_qtl_sum_stats_chr5.txt.gz", package = "cacti"),
#'   out_prefix = out_prefix,
#'   dir_pco = system.file("pco", package = "cacti"),
#'   min_peaks = 1
#' )
#' }
#'
#' @export
cacti_run_chr <- function(
    window_size,
    file_pheno_meta,
    file_pheno,
    file_cov,
    chr,
    qtl_file = NULL,
    file_vcf = NULL,
    file_geno = NULL,
    file_snp_pos = NULL,
    cis_dist = 100000,
    p_threshold = 1.0,
    out_prefix,
    dir_pco = system.file("pco", package = "cacti"),
    min_peaks = 1
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required.")

  # ---- setup paths ----
  out_dir <- dirname(out_prefix)
  if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)

  message("=======================================================")
  message(" CACTI Pipeline")
  message("=======================================================")
  message("Output Directory: ", out_dir)


  tag_window <- gsub("\\s+", "", as.character(window_size))
  file_peak_group  <- paste0(out_prefix, "_peak_group_window", tag_window, ".txt")
  file_peak_group_peaklevel  <- paste0(out_prefix, "_peak_group_window", tag_window, "_peak_as_row.txt")
  file_pheno_cov_residual  <- paste0(out_prefix, "_pheno_cov_residual.txt")
  file_p_peak_group <- file.path(out_dir, paste0(basename(out_prefix), "_pval_window", tag_window, "_", chr, ".txt.gz"))
  matrixqtl_out <- file.path(out_dir, paste0(basename(out_prefix), "_matrixqtl_cis_", chr, ".txt.gz"))

  if (is.null(qtl_file)) {
    if (is.null(file_vcf) && is.null(file_geno)) {
      stop("Provide `qtl_file`, or provide genotype input via `file_vcf` or `file_geno`.")
    }
    message("\n[Step 0/3] Run MatrixEQTL to generate cis QTL summary stats …")
    cacti_matrixqtl_cis(
      file_pheno = file_pheno,
      file_pheno_meta = file_pheno_meta,
      file_cov = file_cov,
      file_qtl_out = matrixqtl_out,
      file_vcf = file_vcf,
      file_geno = file_geno,
      file_snp_pos = file_snp_pos,
      cis_dist = cis_dist,
      p_threshold = p_threshold
    )
    qtl_file <- matrixqtl_out
  }

  message("\n[Step 1/3] Group peaks into non-overlapping windows …")
  cacti_group_peak_window(
    window_size = window_size,
    file_pheno_meta = file_pheno_meta,
    file_peak_group = file_peak_group,
    file_peak_group_peaklevel = file_peak_group_peaklevel
  )

  message("\n[Step 2/3] Residualize phenotype by covariates …")
  cacti_pheno_cov_residual(
    file_pheno = file_pheno,
    file_cov = file_cov,
    file_pheno_cov_residual = file_pheno_cov_residual
  )

  message("\n[Step 3/3] Run pvalue for ", chr, " …")
  cacti_cal_p(
    file_qtl_cis_norm = qtl_file,
    chr = chr,
    file_peak_group = file_peak_group,
    file_pheno_cov_residual = file_pheno_cov_residual,
    file_p_peak_group = file_p_peak_group,
    dir_pco = dir_pco,
    min_peaks = min_peaks
  )

  # ---- final summary message ----
  message("\n=======================================================")
  message("\n CACTI pipeline for ", chr, " completed!\n")
  message("=======================================================")
  message("=======================================================")

  message("----------- Run summary -----------")
  message("Chromosome:", chr)
  message("Window size:", window_size)

  message("\n ----------- Input files -----------")
  message("Peaks BED:", file_pheno_meta)
  message("Phenotype:", file_pheno)
  message("Covariate:", file_cov)
  message("QTL summary stats:", qtl_file)

  message("\n ----------- Output files -----------")
  message("  [1] Grouped peak: ", file_peak_group)
  message("  [2] Grouped peak (peak as row):", file_peak_group_peaklevel)
  message("  [3] Residual phenotype:", file_pheno_cov_residual)
  message("  [4] Pval results (", chr, "): ", file_p_peak_group, "\n")

  message("=======================================================")

  invisible(list(
    file_peak_group = file_peak_group,
    file_peak_group_peaklevel = file_peak_group_peaklevel,
    file_pheno_cov_residual = file_pheno_cov_residual,
    file_p_peak_group = file_p_peak_group
  ))
}
