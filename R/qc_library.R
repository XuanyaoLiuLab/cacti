#' Inspect sequencing/library quality in BAM files
#'
#' An optional, read-only diagnostic. Summarizes reads, mapping and base quality in
#' each supplied BAM without changing the BAM or running a mapping workflow.
#'
#' @param bam_files Named character vector of local BAM file paths. Names are
#'   sample/library IDs and must be unique. An index is not required.
#' @param min_mapq Minimum mapping quality for usable reads (integer 0--254;
#'   default 10). MAPQ 255 or missing means unavailable and never passes.
#' @param exclude_duplicates Exclude reads flagged as duplicates from the
#'   usable-read count? Default `TRUE`; this relies on existing BAM flags.
#' @param duplicate_flags Either `"unknown"` (default) or `"marked"`. Set
#'   `"marked"` only when duplicate marking was performed on all supplied BAMs.
#'   Otherwise the duplicate proportion is `NA`, not evidence of no duplicates.
#' @param min_usable_reads Optional minimum usable-read count for a warning.
#'   Default `NULL`: no study-specific read-depth threshold is imposed.
#' @param warn Emit collected diagnostic warnings? Default `TRUE`.
#' @param yield_size Number of alignment records read per chunk. Default 1e6.
#'
#' @details
#' Counts are reads, not paired-end fragments: each mate is counted separately.
#' `alignment_records` includes all records; all other read counts exclude
#' secondary and supplementary alignments. `total_reads` includes unmapped
#' primary records. `usable_reads` requires mapped primary records, known MAPQ
#' at least `min_mapq`, no QC-fail flag, and (by default) no duplicate flag.
#' These diagnostic criteria do not change CACTI-S counting settings.
#'
#' `mapped_prop` and `usable_prop` use `total_reads` as denominator.
#' `mapq_pass_prop` and `duplicate_prop` use `mapped_reads` as denominator;
#' reads with unavailable MAPQ remain in the mapq-pass denominator but cannot
#' pass. `mean_mapq` excludes unavailable values. Zero denominators give `NA`.
#' `flagged_duplicate_reads` counts mapped primary reads with flag 0x400 even
#' when duplicate-marking provenance is unknown.
#'
#' `mean_base_quality` is the arithmetic mean Phred score across all bases with
#' stored quality, not the mean of per-read means. `prop_bases_q30` is the fraction
#' of those bases with Phred quality at least 30. Both use all primary records,
#' including unmapped, low-MAPQ, duplicate-flagged and QC-failed reads; secondary
#' and supplementary records are excluded. These base-quality summaries are
#' independent of `min_mapq` and `exclude_duplicates`. They do not filter reads.
#' `bases_with_quality` reports their common denominator.
#' `reads_without_base_quality` counts primary records with absent or incomplete
#' stored base qualities, including records with no stored sequence/quality.
#' Missing scores are excluded, not treated as Q0. Both summaries are `NA` when
#' there are no observed quality scores. A mean Phred score is not an arithmetic
#' mean error probability. Base quality describes nucleotide-call confidence;
#' MAPQ describes alignment confidence.
#'
#' Metrics describe only the supplied BAM. Previously removed unmapped,
#' low-quality or duplicate reads cannot be recovered. Even with
#' `duplicate_flags = "marked"`, the duplicate proportion in a filtered BAM
#' is not the original library duplication rate. There is no universal read-depth
#' cutoff across marks or studies. No enrichment or library-complexity metrics
#' are calculated.
#'
#' @return A list with `summary` (one row per BAM), `mapq_distribution`
#'   (counts for MAPQ 0--254 and unavailable 255 among mapped primary reads),
#'   `settings`, `notes`, and `warnings` (character vector). No files are written.
#'   The summary includes `mean_base_quality`, `prop_bases_q30`,
#'   `bases_with_quality`, and `reads_without_base_quality` in addition to the
#'   read-count and mapping-quality metrics described above.
#'
#' @examples
#' bam <- system.file("extdata", "Sample1.bam", package = "cacti")
#' qc <- cacti_qc_library(c(Sample1 = bam), warn = FALSE)
#' qc$summary[, c("sample_id", "total_reads", "mapped_reads", "usable_reads",
#'                 "duplicate_prop")]
#' qc$summary[, c("sample_id", "mean_base_quality", "prop_bases_q30")]
#' qc$warnings
#' @export
cacti_qc_library <- function(
    bam_files, min_mapq = 10, exclude_duplicates = TRUE,
    duplicate_flags = c("unknown", "marked"), min_usable_reads = NULL,
    warn = TRUE, yield_size = 1000000L
) {
  duplicate_flags <- match.arg(duplicate_flags)
  if (!is.character(bam_files) || !length(bam_files) ||
      anyNA(bam_files) || any(!nzchar(bam_files))) {
    stop("bam_files must be a nonempty named character vector of BAM paths.",
         call. = FALSE)
  }
  ids <- names(bam_files)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("bam_files must have unique, nonempty sample/library names.", call. = FALSE)
  }
  if (any(!file.exists(bam_files)) || any(file.info(bam_files)$isdir)) {
    stop("Every bam_files path must identify an existing BAM file.", call. = FALSE)
  }
  valid_number <- function(x, minimum, maximum = Inf, integer = FALSE) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= minimum && x <= maximum && (!integer || x == floor(x))
  }
  if (!valid_number(min_mapq, 0, 254, integer = TRUE)) {
    stop("min_mapq must be an integer from 0 to 254.", call. = FALSE)
  }
  if (!valid_number(yield_size, 1, .Machine$integer.max, integer = TRUE)) {
    stop("yield_size must be a positive integer.", call. = FALSE)
  }
  if (!is.null(min_usable_reads) && !valid_number(min_usable_reads, 0)) {
    stop("min_usable_reads must be NULL or a nonnegative finite number.", call. = FALSE)
  }
  for (value in list(exclude_duplicates, warn)) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("exclude_duplicates and warn must each be TRUE or FALSE.", call. = FALSE)
    }
  }
  ratio <- function(numerator, denominator) {
    if (denominator > 0) numerator / denominator else NA_real_
  }
  scan_one <- function(path, id) {
    bam <- Rsamtools::BamFile(path, yieldSize = as.integer(yield_size))
    open(bam)
    on.exit(close(bam))
    param <- Rsamtools::ScanBamParam(what = c("flag", "mapq", "qual"))
    n_records <- n_total <- n_mapped <- n_usable <- n_dup <- n_qcfail <- 0
    n_pass <- 0
    n_bases_with_quality <- base_quality_sum <- n_bases_q30 <- n_missing_quality <- 0
    hist <- numeric(256L)
    repeat {
      chunk <- Rsamtools::scanBam(bam, param = param)[[1L]]
      if (!length(chunk$flag)) break
      n_records <- n_records + length(chunk$flag)
      primary <- bitwAnd(chunk$flag, 2304L) == 0L # 0x100 | 0x800
      flags <- chunk$flag[primary]
      mapq <- chunk$mapq[primary]
      mapped <- bitwAnd(flags, 4L) == 0L
      duplicate <- bitwAnd(flags, 1024L) != 0L
      qcfail <- bitwAnd(flags, 512L) != 0L
      known_mapq <- !is.na(mapq) & mapq < 255L
      pass_mapq <- known_mapq & mapq >= min_mapq
      usable <- mapped & pass_mapq & !qcfail
      if (exclude_duplicates) usable <- usable & !duplicate
      n_total <- n_total + length(flags)
      n_mapped <- n_mapped + sum(mapped)
      n_usable <- n_usable + sum(usable)
      n_dup <- n_dup + sum(mapped & duplicate)
      n_qcfail <- n_qcfail + sum(mapped & qcfail)
      n_pass <- n_pass + sum(mapped & pass_mapq)
      mq <- mapq[mapped]
      mq[is.na(mq)] <- 255L
      hist <- hist + tabulate(mq + 1L, nbins = 256L)

      # Rsamtools stores Phred scores as ASCII Q+33. Missing BAM quality (255)
      # becomes a space (ASCII 32); absent sequence/quality has zero width.
      # A byte histogram avoids expanding a chunk into per-base integer lists.
      qualities <- Biostrings::BStringSet(chunk$qual[primary])
      quality_hist <- Biostrings::alphabetFrequency(qualities, collapse = TRUE)
      phred_scores <- 0:93
      quality_counts <- as.numeric(quality_hist[phred_scores + 34L])
      n_bases_with_quality <- n_bases_with_quality + sum(quality_counts)
      base_quality_sum <- base_quality_sum + sum(phred_scores * quality_counts)
      n_bases_q30 <- n_bases_q30 + sum(quality_counts[phred_scores >= 30L])
      n_missing_quality <- n_missing_quality + sum(
        Biostrings::width(qualities) == 0L |
          Biostrings::letterFrequency(qualities, " ")[, 1L] > 0L
      )
    }
    known_n <- sum(hist[1:255])
    summary <- data.frame(
      sample_id = id, bam_file = path, alignment_records = n_records,
      total_reads = n_total, mapped_reads = n_mapped,
      mapped_prop = ratio(n_mapped, n_total), usable_reads = n_usable,
      usable_prop = ratio(n_usable, n_total),
      mean_mapq = ratio(sum((0:254) * hist[1:255]), known_n),
      mapq_pass_reads = n_pass, mapq_pass_prop = ratio(n_pass, n_mapped),
      mapq_unavailable_reads = hist[256],
      mean_base_quality = ratio(base_quality_sum, n_bases_with_quality),
      prop_bases_q30 = ratio(n_bases_q30, n_bases_with_quality),
      bases_with_quality = n_bases_with_quality,
      reads_without_base_quality = n_missing_quality,
      flagged_duplicate_reads = n_dup,
      duplicate_prop = if (duplicate_flags == "marked") ratio(n_dup, n_mapped) else NA_real_,
      qc_failed_mapped_reads = n_qcfail, stringsAsFactors = FALSE
    )
    list(summary = summary, histogram = data.frame(
      sample_id = id, mapq = 0:255, available = c(rep(TRUE, 255), FALSE),
      reads = hist, stringsAsFactors = FALSE
    ))
  }
  result <- lapply(seq_along(bam_files), function(i) scan_one(bam_files[[i]], ids[[i]]))
  summary <- do.call(rbind, lapply(result, `[[`, "summary"))
  histogram <- do.call(rbind, lapply(result, `[[`, "histogram"))
  issues <- character()
  if (duplicate_flags == "unknown") {
    issues <- c(issues, paste0("Duplicate-marking status is unknown; duplicate_prop is NA. ",
      "Only existing duplicate flags can be used to exclude reads."))
  }
  for (i in seq_len(nrow(summary))) {
    prefix <- paste0(summary$sample_id[i], ": ")
    if (summary$total_reads[i] == 0) {
      issues <- c(issues, paste0(prefix, "no primary read records."))
    } else if (summary$usable_reads[i] == 0) {
      issues <- c(issues, paste0(prefix, "no reads pass the selected usable-read criteria."))
    }
    if (!is.null(min_usable_reads) && summary$usable_reads[i] < min_usable_reads) {
      issues <- c(issues, paste0(prefix, "usable reads below min_usable_reads (",
                                 min_usable_reads, ")."))
    }
    if (summary$mapq_unavailable_reads[i] > 0) {
      issues <- c(issues, paste0(prefix, "some mapped reads have unavailable MAPQ; ",
                                 "they are excluded from usable reads."))
    }
    if (summary$reads_without_base_quality[i] > 0) {
      issues <- c(issues, paste0(prefix,
        "some primary reads have missing base qualities; summaries use observed scores only."))
    }
  }
  if (warn && length(issues)) warning(paste(issues, collapse = "\n"), call. = FALSE)
  list(summary = summary, mapq_distribution = histogram,
       settings = list(min_mapq = min_mapq, exclude_duplicates = exclude_duplicates,
                       duplicate_flags = duplicate_flags, min_usable_reads = min_usable_reads),
       notes = c("Counts are primary read records, not paired-end fragments.",
                 "Base-quality metrics use observed scores from all primary reads, independent of usable-read filters.",
                 "Metrics describe the supplied BAM only; removed reads cannot be assessed.",
                 "These optional diagnostics do not alter CACTI-S read-counting settings."),
       warnings = issues)
}
