#' Inspect Sample-Level Phenotype Quality
#'
#' Reports missing or non-finite measurements and, for raw counts, total counts
#' for the sample and proportion of zero measurements. Total counts are summed
#' across observed features in the supplied matrix, not all reads in a BAM file.
#' This optional function does not modify data or run automatically during mapping.
#'
#' @param phenotype Numeric R matrix or numeric data frame, with features in rows
#'   and samples in columns. Column names must be sample IDs; row names are
#'   optional. Load tables and remove nonnumeric metadata columns first.
#' @param input_type Either `"raw_counts"` for nonnegative counts (fractional
#'   counts are accepted), or `"normalized"` for normalized/residualized values.
#'   Raw-count diagnostics are not calculated for normalized values.
#' @param min_total_counts Optional nonnegative minimum total counts per sample.
#'   Only applies to raw counts; `NULL` disables this threshold.
#' @param max_zero_prop Optional maximum proportion of zero measurements, between
#'   0 and 1. Only applies to raw counts; `NULL` disables this threshold.
#' @param warn Logical; emit one warning summarizing detected issues. Warnings
#'   remain in the result when `warn = FALSE`.
#'
#' @details `NA` and `NaN` are missing; non-finite values include missing values
#'   and positive/negative infinity. Count totals sum finite observations only;
#'   the zero proportion uses the number of finite observations (`n_observed`)
#'   as its denominator. Both are `NA` when there are no finite observations,
#'   or for normalized inputs. With missing measurements, totals describe only
#'   observed features and should not be interpreted as complete library sizes.
#'   User-specified count thresholds
#'   are ignored for normalized inputs, with a notice in `warnings`. Thresholds
#'   are inspection criteria, not universal quality standards for histone marks.
#'
#' @return A list with `n_samples`, `n_features`, `summary` (one row per sample,
#'   containing `sample_id`, `n_features`, `n_observed`, `n_missing`,
#'   `n_nonfinite`, `prop_missing`, `prop_nonfinite`, `total_counts`, and
#'   `prop_zero`) and `warnings` (a character vector). Missing/non-finite
#'   proportions use all supplied features as their denominator.
#' @export
#' @examples
#' y <- matrix(c(0, 4, 2, 6, 3, 8), nrow = 2,
#'             dimnames = list(c("peak1", "peak2"), c("s1", "s2", "s3")))
#' cacti_qc_samples(y, input_type = "raw_counts")$summary
cacti_qc_samples <- function(phenotype,
                             input_type = c("raw_counts", "normalized"),
                             min_total_counts = NULL, max_zero_prop = NULL,
                             warn = TRUE) {
  .cacti_qc_check_warn(warn)
  input_type <- match.arg(input_type)
  x <- .cacti_qc_matrix(phenotype, "phenotype")
  .cacti_qc_threshold(min_total_counts, "min_total_counts", nullable = TRUE)
  .cacti_qc_threshold(max_zero_prop, "max_zero_prop", upper = 1, nullable = TRUE)
  finite <- is.finite(x)
  if (input_type == "raw_counts") .cacti_qc_nonnegative(x, finite)
  n_missing <- colSums(is.na(x))
  n_nonfinite <- colSums(!finite)
  n_observed <- colSums(finite)
  totals <- zeros <- rep(NA_real_, ncol(x))
  issues <- character()
  if (any(n_nonfinite > 0)) {
    issues <- c(issues, sprintf("%d sample(s) contain missing or non-finite measurements.",
                                sum(n_nonfinite > 0)))
  }
  if (input_type == "raw_counts") {
    observed <- x
    observed[!finite] <- NA_real_
    totals <- colSums(observed, na.rm = TRUE)
    zeros <- colSums(observed == 0, na.rm = TRUE) / n_observed
    totals[n_observed == 0] <- NA_real_
    zeros[n_observed == 0] <- NA_real_
    overflow <- n_observed > 0 & !is.finite(totals)
    if (any(overflow)) {
      totals[overflow] <- NA_real_
      issues <- c(issues, "Some raw-count totals exceeded the numeric range and were set to NA.")
    }
    if (any(zeros == 1, na.rm = TRUE)) {
      issues <- c(issues, sprintf("%d sample(s) have zero counts for every observed feature.",
                                  sum(zeros == 1, na.rm = TRUE)))
    }
    if (!is.null(min_total_counts) && any(totals < min_total_counts, na.rm = TRUE)) {
      issues <- c(issues, sprintf("%d sample(s) have total counts below min_total_counts.",
                                  sum(totals < min_total_counts, na.rm = TRUE)))
    }
    if (!is.null(max_zero_prop) && any(zeros > max_zero_prop, na.rm = TRUE)) {
      issues <- c(issues, sprintf("%d sample(s) have a zero proportion above max_zero_prop.",
                                  sum(zeros > max_zero_prop, na.rm = TRUE)))
    }
  } else if (!is.null(min_total_counts) || !is.null(max_zero_prop)) {
    issues <- c(issues, "Raw-count thresholds are not applied to normalized inputs.")
  }
  summary <- data.frame(
    sample_id = colnames(x), n_features = nrow(x), n_observed = n_observed,
    n_missing = n_missing,
    n_nonfinite = n_nonfinite, prop_missing = n_missing / nrow(x),
    prop_nonfinite = n_nonfinite / nrow(x), total_counts = totals,
    prop_zero = zeros, stringsAsFactors = FALSE, row.names = NULL
  )
  .cacti_qc_emit(issues, warn)
  list(n_samples = ncol(x), n_features = nrow(x), summary = summary, warnings = issues)
}

#' Inspect Feature-Level Phenotype Quality
#'
#' Reports missingness, variation, and (for raw counts) abundance and detection.
#' Returns a descriptive table, not pass/fail classifications. Optional warning
#' criteria do not remove features. Independent of CACTI mapping functions.
#'
#' @param phenotype Numeric R matrix or numeric data frame, with features in rows
#'   and samples in columns. Row names must be feature IDs and column names must
#'   be sample IDs. Load tables and remove nonnumeric metadata columns first.
#' @param input_type Either `"raw_counts"` for nonnegative counts (fractional
#'   counts are accepted), or `"normalized"` for normalized/residualized values.
#' @param min_count Nonnegative count threshold defining detection: a finite
#'   raw count is detected when it is greater than or equal to `min_count`.
#'   Ignored for normalized inputs.
#' @param min_prop Warn when the detected proportion of observed samples is
#'   below this value (between 0 and 1). Ignored for normalized inputs.
#' @param min_variance Warn when variance is not greater than this nonnegative
#'   value. The default warns about constant features; nothing is removed.
#' @param max_missing_prop Warn when the proportion of non-finite measurements
#'   exceeds this value (between 0 and 1); includes `NA`, `NaN`, and infinity.
#' @param warn Logical; emit one warning summarizing detected issues. Warnings
#'   remain in the result when `warn = FALSE`.
#'
#' @details Variance and raw-count mean abundance use finite observations only.
#'   Variance is `NA` with fewer than two finite observations. Detection and zero
#'   proportions use the number of finite observations (`n_observed`) as their
#'   denominator. Count-based metrics are `NA` if nothing is observed, or for
#'   normalized inputs. Missing/non-finite proportions use all supplied samples
#'   as their denominator. Warning thresholds are inspection criteria, not
#'   universal quality cutoffs or evidence of calibration or statistical power.
#'
#' @return A list containing:
#' \itemize{
#'   \item `n_samples` and `n_features`: dimensions of the supplied matrix.
#'   \item `summary`: one row per feature, containing `feature_id`, `n_samples`,
#'     `n_observed`, `n_missing`, `n_nonfinite`, `prop_missing`, `prop_nonfinite`,
#'     `variance`, `mean_count`, `n_detected`, `prop_detected`, and `prop_zero`.
#'   \item `warnings`: character vector describing detected issues.
#' }
#' @export
#' @examples
#' y <- matrix(c(0, 4, 2, 6, 3, 8), nrow = 2,
#'             dimnames = list(c("peak1", "peak2"), c("s1", "s2", "s3")))
#' qc <- cacti_qc_features(y, input_type = "raw_counts", min_count = 1,
#'                         min_prop = 0.5)
#' qc$summary
#' qc$n_features
cacti_qc_features <- function(phenotype,
                              input_type = c("raw_counts", "normalized"),
                              min_count = 1, min_prop = 0.2, min_variance = 0,
                              max_missing_prop = 0, warn = TRUE) {
  .cacti_qc_check_warn(warn)
  input_type <- match.arg(input_type)
  x <- .cacti_qc_matrix(phenotype, "phenotype", require_rows = TRUE)
  .cacti_qc_threshold(min_count, "min_count")
  .cacti_qc_threshold(min_prop, "min_prop", upper = 1)
  .cacti_qc_threshold(min_variance, "min_variance")
  .cacti_qc_threshold(max_missing_prop, "max_missing_prop", upper = 1)
  finite <- is.finite(x)
  if (input_type == "raw_counts") .cacti_qc_nonnegative(x, finite)
  n_missing <- rowSums(is.na(x))
  n_nonfinite <- rowSums(!finite)
  n_observed <- rowSums(finite)
  variances <- apply(x, 1L, function(values) {
    observed <- values[is.finite(values)]
    if (length(observed) >= 2L) stats::var(observed) else NA_real_
  })
  means <- detected <- prop_detected <- zeros <- rep(NA_real_, nrow(x))
  if (input_type == "raw_counts") {
    observed <- x
    observed[!finite] <- NA_real_
    means <- rowMeans(observed, na.rm = TRUE)
    means[!is.finite(means)] <- NA_real_
    detected <- rowSums(finite & x >= min_count, na.rm = TRUE)
    detected[n_observed == 0] <- NA_real_
    prop_detected <- detected / n_observed
    zeros <- rowSums(observed == 0, na.rm = TRUE) / n_observed
    zeros[n_observed == 0] <- NA_real_
  }
  summary <- data.frame(
    feature_id = rownames(x), n_samples = ncol(x), n_observed = n_observed,
    n_missing = n_missing,
    n_nonfinite = n_nonfinite, prop_missing = n_missing / ncol(x),
    prop_nonfinite = n_nonfinite / ncol(x), variance = variances,
    mean_count = means, n_detected = detected, prop_detected = prop_detected,
    prop_zero = zeros,
    stringsAsFactors = FALSE, row.names = NULL
  )
  issues <- character()
  if (any(n_nonfinite / ncol(x) > max_missing_prop)) {
    issues <- c(issues, sprintf("%d feature(s) have a non-finite proportion above max_missing_prop.",
                                sum(n_nonfinite / ncol(x) > max_missing_prop)))
  }
  low_variance <- !is.finite(variances) | variances <= min_variance
  if (any(low_variance)) {
    issues <- c(issues, sprintf(
      "%d feature(s) have insufficient finite observations or variance not above min_variance.",
      sum(low_variance)))
  }
  if (input_type == "raw_counts" && any(prop_detected < min_prop, na.rm = TRUE)) {
    issues <- c(issues, sprintf("%d feature(s) have an observed detection proportion below min_prop.",
                                sum(prop_detected < min_prop, na.rm = TRUE)))
  }
  .cacti_qc_emit(issues, warn)
  list(n_samples = ncol(x), n_features = nrow(x), summary = summary, warnings = issues)
}

# Shared validation is deliberately limited to the documented matrix interface.
.cacti_qc_matrix <- function(x, name, require_rows = FALSE) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, function(column) is.numeric(column) && !is.complex(column),
                    logical(1)))) {
      stop(name, " must contain only numeric columns; remove metadata columns first.",
           call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) == 0L || ncol(x) == 0L) {
    stop(name, " must be a nonempty numeric matrix or numeric data frame.", call. = FALSE)
  }
  valid_ids <- function(ids) !is.null(ids) && !anyNA(ids) && all(nzchar(trimws(ids)))
  if (!valid_ids(colnames(x))) {
    stop(name, " must have nonempty sample IDs as column names.", call. = FALSE)
  }
  if (require_rows && !valid_ids(rownames(x))) {
    stop(name, " must have nonempty feature IDs as row names.",
         call. = FALSE)
  }
  x
}

.cacti_qc_check_warn <- function(warn) {
  if (!is.logical(warn) || length(warn) != 1L || is.na(warn)) {
    stop("warn must be TRUE or FALSE.", call. = FALSE)
  }
}

.cacti_qc_threshold <- function(value, name, upper = Inf, nullable = FALSE) {
  if (nullable && is.null(value)) return(invisible(NULL))
  if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
      !is.finite(value) || value < 0 || value > upper) {
    stop(name, " must be a finite number between 0 and ", upper,
         if (nullable) ", or NULL" else "", ".", call. = FALSE)
  }
  invisible(NULL)
}

.cacti_qc_nonnegative <- function(x, finite) {
  if (any(x[finite] < 0)) {
    stop("raw_counts must not contain negative finite values; use input_type = 'normalized' for normalized data.",
         call. = FALSE)
  }
}

.cacti_qc_emit <- function(issues, warn) {
  if (warn && length(issues)) warning(paste(issues, collapse = "\n"), call. = FALSE)
  invisible(NULL)
}
