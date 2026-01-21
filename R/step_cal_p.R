#' Run CACTI association test for a single chromosome
#'
#' @param file_qtl_cis_norm Path to cis QTL summary stats.
#' Must contain columns: phe_id, var_id, z.
#' @param chr Chromosome label (e.g., "chr21"). Used to subset peak groups.
#' @param file_peak_group Path to peak-window metadata (groups) file. Produced by step_group_peak_window.
#' @param file_pheno_cov_residual Path to residualized phenotype matrix file. Produced by step_pheno_cov_residual.
#' @param file_p_peak_group Output path for per-window pvalues.
#' Columns: group, snp, pval.
#' @param dir_pco Directory containing association test helpers: ModifiedPCOMerged_acat.R, liu.R, liumod.R, davies.R, qfc.so.
#' @param min_peaks Minimum number of peaks required in a group to run the multivariate PCO test (>= min_peaks -> PCO; < min_peaks -> univariate p).
#' @return Invisibly returns the pvalue tibble; writes `file_p_peak_group`.
#' @export
cacti_cal_p <- function(
    file_qtl_cis_norm, chr,
    file_peak_group, file_pheno_cov_residual,
    file_p_peak_group,
    dir_pco = system.file("pco", package = "cacti"),
    min_peaks = 2
    ) {
  suppressPackageStartupMessages(requireNamespace("data.table"))
  suppressPackageStartupMessages(requireNamespace("dplyr"))
  suppressPackageStartupMessages(requireNamespace("tidyr"))
  suppressPackageStartupMessages(requireNamespace("tibble"))
  suppressPackageStartupMessages(requireNamespace("stringr"))
  suppressPackageStartupMessages(requireNamespace("ACAT"))

  ## -- Load PCO helpers --
  if (!is.null(dir_pco)) {
    srcs <- file.path(dir_pco, c("ModifiedPCOMerged_acat.R","liu.R","liumod.R","davies.R"))
    for (s in srcs) if (file.exists(s)) source(s)
    so_path <- file.path(dir_pco, "qfc.so")
    if (file.exists(so_path)) try(dyn.load(so_path), silent = TRUE)
  }
  if (!exists("ModifiedPCOMerged_acat")) {
    stop("ModifiedPCOMerged_acat not found. Provide 'dir_pco' with helper scripts or attach the function in your session.")
  }

  ## -- Read inputs --
  qtl_cis_norm <- data.table::fread(file_qtl_cis_norm, header = TRUE)
  if(!all(c('phe_id', 'var_id', 'z') %in% colnames(qtl_cis_norm))) stop("cis QTL summary stats must contain columns: 'phe_id', 'var_id', 'z'.")

  peak_group <- data.table::fread(file_peak_group, header = TRUE)

  pheno_cov_residual <- data.table::fread(file_pheno_cov_residual, header = TRUE)
  pheno_cov_residual_mat <- dplyr::select(pheno_cov_residual, -ID) |> as.matrix() |> t()
  colnames(pheno_cov_residual_mat) <- pheno_cov_residual$ID


  ## -- Build peak group indices --
  df_peak_group <- dplyr::filter(peak_group, n_pid_group >= !!min_peaks) |>
    dplyr::select(group, group_chr, n_pid_group, pid_group) |>
    tidyr::separate_rows(pid_group, sep = ";") |>
    dplyr::arrange(group_chr, n_pid_group)
  group_list <- dplyr::filter(df_peak_group, group_chr == !!chr) |>
    dplyr::distinct(group) |>
    dplyr::pull(group)
  n_group <- length(group_list)


  ## -- Single-peak groups processed with univariate p-values --
  df_peak_group_single <- dplyr::filter(peak_group, n_pid_group < !!min_peaks) |>
    dplyr::select(group, group_chr, n_pid_group, pid_group) |>
    tidyr::separate_rows(pid_group, sep = ";") |>
    dplyr::arrange(group_chr, n_pid_group)
  group_list_single <- dplyr::filter(df_peak_group_single, group_chr == !!chr) |>
    dplyr::distinct(group) |>
    dplyr::pull(group)
  n_peak_group_single <- length(group_list_single)


  ## -- Main loop over groups --
  p_peak_group_list <- vector("list", length = n_group + n_peak_group_single)
  names(p_peak_group_list) <- c(group_list, group_list_single)
  k <- 0

  # Multi-peak groups: PCO
  for (group in group_list) {
    k <- k + 1
    if (k %% 100 == 0) cat(k, "-th group is running, out of", n_group + n_peak_group_single, "group, on", chr, ".\n")

    tmp_phe_id <- dplyr::filter(df_peak_group, group == !!group) |>
      dplyr::pull(pid_group)

    z_mat <- dplyr::filter(qtl_cis_norm, phe_id %in% !!tmp_phe_id) |>
      dplyr::select(phe_id, var_id, z) |>
      tidyr::pivot_wider(names_from = phe_id, values_from = z) |>
      tidyr::drop_na() |>
      tibble::column_to_rownames(var = "var_id") |>
      as.matrix()

    if (nrow(z_mat) == 0 || ncol(z_mat) < min_peaks) {
      p_peak_group_list[[group]] <- NULL
      next
    }

    Sigma <- stats::cor(pheno_cov_residual_mat[, colnames(z_mat), drop = FALSE])

    vals <- ModifiedPCOMerged_acat(Z.mat = z_mat, Sigma = Sigma)

    p_peak_group_list[[group]] <- tibble::enframe(vals, name = "snp", value = "pval") |>
      dplyr::mutate(group = !!group)
  }



  # Single-peak groups: use univariate p from z
  for (group in group_list_single) {
    k <- k + 1
    if (k %% 100 == 0) cat(k, "-th group is running, out of", n_group + n_group, "group, on", chr, ".\n")

    tmp_phe_id <- dplyr::filter(df_peak_group_single, group == !!group) |>
      dplyr::pull(pid_group)

    p_peak_group_list[[group]] <- dplyr::filter(qtl_cis_norm, phe_id %in% !!tmp_phe_id) |>
      dplyr::select(snp = var_id, z) |>
      tidyr::drop_na() |>
      dplyr::mutate(
        pval = 1 - stats::pchisq(z^2, df = 1),
        z = NULL,
        group = !!group
      )
  }

  out <- dplyr::bind_rows(p_peak_group_list) |>
    dplyr::relocate(group, .before = dplyr::everything())

  data.table::fwrite(out, file = file_p_peak_group, quote = FALSE, sep = "\t", col.names = TRUE)
  invisible(out)
}
