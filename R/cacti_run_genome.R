#' Run CACTI peak-window pipeline genome-wide and add FDR
#'
#' This function is a convenience wrapper that runs the CACTI peak-window pipeline
#' across multiple chromosomes and then computes window-level FDR across all
#' peak windows and chromosomes.
#'
#' For each chromosome in `chrs`, it calls [cacti_run_chr()] and
#' collects the per-window p-value files. It then calls [cacti_add_fdr()]
#' once, aggregating all chromosomes to obtain q-values for the top-hit p-values
#' in each window.
#'
#' @inheritParams cacti_run_chr
#' @param chrs Character vector of chromosome labels (e.g., `paste0("chr", 1:22)`).
#' @param qtl_files Either:
#'   \itemize{
#'     \item a character vector of the same length as `chrs`, giving the cis-QTL
#'           file path for each chromosome in order; or
#'     \item a single template string containing the placeholder `"{chr}"`,
#'           e.g., `"extdata/test_qtl_sum_stats_{chr}.txt.gz"`. In that case,
#'           the placeholder is replaced by each element of `chrs`.
#'   }
#' @param file_fdr_out Optional output path for the FDR-added window-level file.
#'   If `NULL`, a default filename is constructed from `out_prefix` and `window_size`.
#'
#' @return Invisibly returns a named list of output paths with elements:
#'   \describe{
#'     \item{file_peak_group}{Path to the window-level group.}
#'     \item{file_peak_group_peaklevel}{Path to the peak-level group.}
#'     \item{file_pheno_residual}{Path to the residualized phenotype matrix.}
#'     \item{file_p_peak_group}{Path to the per-window p-value file for all chromosome.}
#'     \item{file_fdr_out}{Path to the FDR-added window-level result file.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example (chr5 as all chr)
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
#' qtl_file <- system.file(
#'   "extdata", "test_qtl_sum_stats_chr5.txt.gz",
#'   package = "cacti"
#' )
#'
#' out_prefix <- tempfile("cacti_genome_")
#'
#' res <- cacti_run_genome(
#'   window_size = "50kb",
#'   file_pheno_meta = file_pheno_meta,
#'   file_pheno = file_pheno,
#'   file_cov = file_cov,
#'   chrs = "chr5",
#'   qtl_files = qtl_file,
#'   out_prefix = out_prefix,
#'   dir_pco = system.file("pco", package = "cacti"),
#'   min_peaks = 2,
#'   file_fdr_out = file.path(tempdir(), "cacti_fdr_chr5.txt.gz")
#' )
#' }
#'
#' @export
cacti_run_genome <- function(
    window_size,
    file_pheno_meta,
    file_pheno,
    file_cov,
    chrs,
    qtl_files,
    out_prefix,
    dir_pco = system.file("pco", package = "cacti"),
    min_peaks = 2,
    file_fdr_out   = NULL
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required.")

  # resolve qtl_files
  if (length(qtl_files) == 1L && grepl("\\{chr\\}", qtl_files)) {
    qtl_files <- vapply(
      chrs,
      function(cc) gsub("\\{chr\\}", cc, qtl_files),
      FUN.VALUE = character(1L)
    )
  }
  if (length(qtl_files) != length(chrs)) {
    stop("`qtl_files` must be either length 1 with '{chr}' placeholder, ",
         "or the same length as `chrs`.")
  }

  # setup paths
  out_dir <- dirname(out_prefix)
  if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)

  message("=======================================================")
  message(" CACTI Pipeline Genome-Wide")
  message("=======================================================")
  message("Output Directory: ", out_dir)

  message("--------------------------------------------------------")
  message("=== [I/III] Running CACTI peak-window pipeline genome-wide ===")
  message("--------------------------------------------------------")

  tag_window <- gsub("\\s+", "", as.character(window_size))
  file_peak_group  <- paste0(out_prefix, "_peak_group_window", tag_window, ".txt")
  file_peak_group_peaklevel  <- paste0(out_prefix, "_peak_group_window", tag_window, "_peak_as_row.txt")
  file_pheno_cov_residual  <- paste0(out_prefix, "_pheno_cov_residual.txt")


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

  pval_files <- character(length(chrs))
  for (i in seq_along(chrs)) {
    this_chr <- chrs[i]
    this_qtl <- qtl_files[i]

    message("\n[Step 3/3] Run pvalue for ", this_chr, " …")

    file_p_peak_group <- file.path(out_dir, paste0(basename(out_prefix), "_pval_window", tag_window, "_", this_chr, ".txt.gz"))
    cacti_cal_p(
      file_qtl_cis_norm = this_qtl,
      chr = this_chr,
      file_peak_group = file_peak_group,
      file_pheno_cov_residual = file_pheno_cov_residual,
      file_p_peak_group = file_p_peak_group,
      dir_pco = dir_pco,
      min_peaks = min_peaks
    )
    pval_files[i] <- file_p_peak_group


    # ---- final summary message ----
    message("\n pipeline for ", this_chr, " completed!\n")

    message("\n ----------- Input files -----------")
    message("Peaks BED:", file_pheno_meta)
    message("Phenotype:", file_pheno)
    message("Covariate:", file_cov)
    message("QTL summary stats:", this_qtl)

    message("\n ----------- Output files -----------")
    message("  [1] Grouped peak: ", file_peak_group)
    message("  [2] Grouped peak (peak as row):", file_peak_group_peaklevel)
    message("  [3] Residual phenotype:", file_pheno_cov_residual)
    message("  [4] Pval results (", this_chr, "): ", file_p_peak_group, "\n")

  }

  if (is.null(file_fdr_out)) {
    tag_window <- gsub("\\s+", "", as.character(window_size))
    file_fdr_out <- file.path(
      dirname(out_prefix),
      paste0(basename(out_prefix), "_pval_window", tag_window, "_fdr_added.txt.gz")
    )
  }

  message("--------------------------------------------------------")
  message("\n=== [II/III] Adding FDR across windows and chromosomes ===")
  message("--------------------------------------------------------")
  cacti_add_fdr(
    file_all_pval = pval_files,
    file_fdr_out = file_fdr_out
  )


  message("--------------------------------------------------------")
  message("\n=== [III/III] Genome-wide CACTI peak-window pipeline completed with FDR! ===\n")
  message("--------------------------------------------------------")

  message("----------- Run summary -----------")
  message("  [1] Pvalue output: ", paste(pval_files, collapse = ";"))
  message("  [2] FDR output: ", file_fdr_out)

  invisible(list(
    file_peak_group = file_peak_group,
    file_peak_group_peaklevel = file_peak_group_peaklevel,
    file_pheno_cov_residual = file_pheno_cov_residual,
    file_p_peak_group = pval_files,
    file_fdr_out = file_fdr_out
  ))
}
