#' Group adjacent peaks by non-overlapping windows
#'
#' Group peaks into non-overlapping genomic windows of size `window_size`,
#' using peak coordinates from a BED-like metadata table.
#'
#' @param window_size Character like `"50kb"` (kb units) or numeric (bp).
#' @param file_pheno_meta Path to input phenotype meta table (with header). Columns are -
#' - `phe_chr`: chromosome, with "chr" prefix (e.g. "chr1")
#' - `phe_from`: peak start (integer, < `phe_to`)
#' - `phe_to`: peak end (integer)
#' - `phe_id`: peak identifier (unique per row)
#' @param file_peak_group Output path for window-level table. Columns are -
#' group, group_chr, group_from, group_to, pid_group, n_pid_group
#' @param file_peak_group_peaklevel Output path for peak-level table. Columns are -
#' phe_id, phe_chr, phe_from, phe_to, group, n_pid_group
#'
#' @return Invisibly returns the grouped window-level table; writes `file_peak_group`.
#' @export
cacti_group_peak_window <- function(window_size, file_pheno_meta, file_peak_group, file_peak_group_peaklevel) {
  # packages namespace
  suppressPackageStartupMessages(requireNamespace("data.table"))
  suppressPackageStartupMessages(requireNamespace("dplyr"))
  suppressPackageStartupMessages(requireNamespace("tidyr"))
  suppressPackageStartupMessages(requireNamespace("tibble"))
  suppressPackageStartupMessages(requireNamespace("stringr"))

  # window size
  if (is.character(window_size) && stringr::str_detect(window_size, "^\\d+[.]?\\d*[kK][bB]$")) {
    window_size <- as.numeric(stringr::str_extract(window_size, "^\\d+[.]?\\d*")) * 1000
  } else if (!is.numeric(window_size)) {
    stop("Specify window size in kb like '50kb' or provide a numeric bp value.")
  }

  # read pheno meta
  pheno_meta <- data.table::fread(file_pheno_meta, header = TRUE)
  colnames(pheno_meta) <- c('phe_id', 'phe_chr', 'phe_from', 'phe_to')

  # add chr prefix if missing
  if (!all(startsWith(pheno_meta$phe_chr, "chr")))
    pheno_meta$phe_chr <- ifelse(startsWith(pheno_meta$phe_chr, "chr"), pheno_meta$phe_chr, paste0("chr", pheno_meta$phe_chr))

  # autosomes only
  pheno_meta <- dplyr::filter(pheno_meta, stringr::str_detect(phe_chr, '^chr\\d+$'))

  # parse coordinates
  if (any(pheno_meta$phe_from >= pheno_meta$phe_to)) stop("feature start must be < end.")

  # group by non-overlapping windows
  df_peak_group <- pheno_meta |>
    dplyr::group_by(phe_chr) |>
    dplyr::mutate(
      group = cut(
        phe_from,
        breaks = c(seq(0, max(phe_from)+1, by = window_size), Inf),
        labels = paste0(unique(phe_chr), "_G", seq_along(seq(0, max(phe_from)+1, by = window_size))),
        right = FALSE
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      group_chr = dplyr::first(phe_chr),
      group_from = min(phe_from),
      group_to = max(phe_to),
      pid_group = paste(phe_id, collapse = ";"),
      n_pid_group = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(group_chr, group_from) |>
    dplyr::mutate(group = NULL) |>
    tibble::rowid_to_column(var = "group") |>
    dplyr::mutate(group = paste0("G", group))

  # peak-level of grouping
  peak_meta_w_group <- tidyr::separate_rows(df_peak_group, pid_group, sep = ";") |>
    dplyr::select(group, phe_id = pid_group, n_pid_group) |>
    dplyr::right_join(pheno_meta, by = c('phe_id')) |>
    dplyr::select(phe_id, phe_chr, phe_from, phe_to, group, n_pid_group)

  # write out
  data.table::fwrite(df_peak_group, file = file_peak_group, sep = '\t', quote = FALSE, col.names = TRUE)
  data.table::fwrite(peak_meta_w_group, file = file_peak_group_peaklevel, sep = '\t', quote = FALSE, col.names = TRUE)

  invisible(peak_meta_w_group)
}
