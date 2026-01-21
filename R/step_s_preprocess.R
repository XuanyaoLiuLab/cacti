#' Preprocess CACTI-S Segments
#'
#' Performs QC, filtering, and normalization
#'
#' **Processing Steps:**
#' \enumerate{
#'   \item \strong{QC:} Remove segments with 0 reads in all samples.
#'   \item \strong{CPM:} Transform counts to log2(CPM + 0.5).
#'   \item \strong{Filter:} Subset segments based on \code{filter_mode}.
#'   \item \strong{Standardize:} Z-score normalization per feature (Mean=0, SD=1).
#'   \item \strong{RankNorm:} Inverse Normal Transform per sample.
#' }
#'
#' **Filtering Logic:**
#' \itemize{
#'   \item \strong{cpm_threshold} (Default): Keeps segments with peak intensity > \code{min_cpm}
#'     in > \code{min_prop} of samples. Standard method to remove noise.
#'   \item \strong{match_overlap_count}: Determines $N$, the number of segments overlapping
#'     regions in \code{filter_bed}. Then keeps the top $N$ segments by average intensity.
#'   \item \strong{prop}: Keeps the top \code{prop_top} proportion of segments by average intensity.
#' }
#'
#' @param file_raw_counts Path to featureCounts output (from \code{cacti_s_count}).
#' @param file_out Output path for the processed phenotype file.
#' @param file_out_meta Output path for the meta file of processed phenotype.
#' @param filter_mode Strategy to select segments: "cpm_threshold" (default), "match_overlap_count", or "prop".
#' @param filter_bed (Required for "match_overlap_count") Path to a BED file (e.g. peaks or QTL regions).
#' @param prop_top (For "prop") Proportion of segments to keep (0.0 to 1.0). Default 0.2.
#' @param min_cpm (For "cpm_threshold") Minimum CPM threshold (default 1).
#' @param min_prop (For "cpm_threshold") Minimum proportion of samples (default 0.2).
#'
#' @return Invisibly returns the processed data frame. Writes \code{file_out} to disk.
#' @export
#' @examples
#' \dontrun{
#' # Locate the bundled test data
#' filecounts <- system.file("extdata/test_results", "test_raw_counts.txt", package = "cacti")
#' filepeaks  <- system.file("extdata", "test_peaks.bed", package = "cacti")
#'
#' # Define output file
#' fileout <- "inst/extdata/test_results/test_pheno_norm.txt"
#' fileout <- "inst/extdata/test_results/test_pheno_norm_meta.txt"
#'
#' # Check if files exist (only true if package is installed)
#' if (nchar(filecounts) > 0) {
#'
#'   # --- Example 1: Default Filtering (CPM Threshold) ---
#'   cacti_s_preprocess(
#'     file_raw_counts = filecounts,
#'     file_out = fileout,
#'     file_out_meta = fileout_meta,
#'     filter_mode = "cpm_threshold",
#'     min_cpm = 10,
#'     min_prop = 0.5
#'   )
#'   head(read.table(fileout, header = TRUE))
#'
#'   # --- Example 2: Match Overlap Count ---
#'   # Logic: If 2 segments overlap peaks, keep the top 2 segments by intensity.
#'   cacti_s_preprocess(
#'     file_raw_counts = filecounts,
#'     file_out = fileout,
#'     file_out_meta = fileout_meta,
#'     filter_mode = "match_overlap_count",
#'     filter_bed = filepeaks
#'   )
#'   head(read.table(fileout, header = TRUE))
#' }
#' }
cacti_s_preprocess <- function(
    file_raw_counts,
    file_out,
    file_out_meta,
    filter_mode = c("cpm_threshold", "match_overlap_count", "prop"),
    filter_bed = NULL,
    prop_top = 0.2,
    min_cpm = 1,
    min_prop = 0.2
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required.")
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR required.")

  # Default is now "cpm_threshold" (the first element)
  filter_mode <- match.arg(filter_mode)

  message("=== Starting Preprocessing ===")
  message("  Input: ", file_raw_counts)
  message("  Filter Mode: ", filter_mode)

  # --- 1. Read Data & Basic QC ---
  raw <- data.table::fread(file_raw_counts)

  # Check standard featureCounts columns
  required_cols <- c("GeneID", "Chr", "Start", "End", "Strand", "Length")
  if (!all(required_cols %in% colnames(raw))) {
    stop("Input file does not look like featureCounts output. Missing columns: ",
         paste(setdiff(required_cols, colnames(raw)), collapse=", "))
  }


  # --- Clean Sample Names ---
  # featureCounts outputs columns like "/path/to/Sample1.bam"
  # We want just "Sample1" to match VCF
  current_cols <- colnames(raw)
  sample_cols_idx <- which(!current_cols %in% required_cols)
  if (length(sample_cols_idx) == 0) stop("No sample columns found in input file.")

  # 1. Remove directory paths (basename)
  # 2. Remove extensions (.bam, .sam, .cram, .txt)
  cleaned_names <- basename(current_cols[sample_cols_idx])
  cleaned_names <- sub("\\.(bam|sam|cram|txt)$", "", cleaned_names, ignore.case = TRUE)
  colnames(raw)[sample_cols_idx] <- cleaned_names

  sample_cols <- cleaned_names
  message("  Cleaned sample names: ", paste(head(sample_cols, 3), collapse=", "), "...")

  counts_mat <- as.matrix(raw[, sample_cols, with = FALSE])
  rownames(counts_mat) <- raw$GeneID

  message("  Initial Features: ", nrow(counts_mat))
  message("  Initial Samples:  ", ncol(counts_mat))


  # [Step 1: QC] Remove segments with 0 reads across ALL samples
  keep_nonzero <- rowSums(counts_mat) > 0
  counts_mat <- counts_mat[keep_nonzero, ]
  raw_subset <- raw[keep_nonzero, ]

  message("  Features after removing all-zeros: ", nrow(counts_mat))

  # --- 2. Normalization (CPM) ---
  # log2(CPM + 0.5)
  cpm_mat <- edgeR::cpm(counts_mat, log = TRUE, prior.count = 0.5)

  # Calculate Average Intensity
  ave_cpm <- rowMeans(cpm_mat)

  # --- 3. Filtering ---
  ids_to_keep <- rownames(cpm_mat)

  if (filter_mode == "match_overlap_count") {
    # Logic:
    # 1. Count N segments that overlap with filter_bed.
    # 2. Keep Top N segments by Intensity.

    if (is.null(filter_bed)) stop("Argument 'filter_bed' is required for 'match_overlap_count' mode.")
    message("  Calculating overlap with: ", filter_bed)

    # Load BED
    peaks <- data.table::fread(filter_bed, select = 1:3)
    colnames(peaks) <- c("chr", "start", "end")

    # Prepare Segments
    segs <- raw_subset[, c("GeneID", "Chr", "Start", "End"), with = FALSE]

    # Find overlaps
    data.table::setkey(peaks, chr, start, end)
    data.table::setkey(segs, Chr, Start, End)
    overlaps <- data.table::foverlaps(segs, peaks, type = "any", nomatch = 0L)
    overlapping_ids <- unique(overlaps$GeneID)

    n_overlap <- length(overlapping_ids)
    message("  Segments physically overlapping peaks: ", n_overlap)

    if (n_overlap == 0) stop("No segments overlapped with the provided BED file.")

    # Select Top N
    message("  Selecting top ", n_overlap, " segments by average intensity...")
    ids_to_keep <- names(sort(ave_cpm, decreasing = TRUE))[1:n_overlap]

  } else if (filter_mode == "prop") {
    # Logic: Keep top X% by Intensity
    if (prop_top <= 0 || prop_top > 1) stop("prop_top must be between 0 and 1.")

    n_keep <- floor(nrow(cpm_mat) * prop_top)
    message("  Filtering top ", prop_top * 100, "% (", n_keep, ") features by average CPM.")

    ids_to_keep <- names(sort(ave_cpm, decreasing = TRUE))[1:n_keep]

  } else if (filter_mode == "cpm_threshold") {
    # Logic: CPM thresholding
    message("  Filtering by intensity: CPM > ", min_cpm, " in > ", min_prop*100, "% samples")

    raw_cpm <- edgeR::cpm(counts_mat, log = FALSE)
    pass_freq <- rowMeans(raw_cpm > min_cpm)
    ids_to_keep <- rownames(cpm_mat)[pass_freq > min_prop]
  }

  # Apply Filter
  message("  Features kept: ", length(ids_to_keep))
  if (length(ids_to_keep) == 0) stop("Filtering removed all features!")

  mat_filtered <- cpm_mat[ids_to_keep, ]

  # --- 4. Standardization (Across Samples) ---
  # Transform each feature to Mean=0, SD=1
  message("  Standardizing (Z-score per feature)...")
  mat_std <- t(scale(t(mat_filtered)))

  # --- 5. Rank Normalization (Across Features) ---
  # Transform each sample to Standard Normal distribution
  message("  Rank Normalizing (INT per sample)...")

  rank_norm <- function(u) {
    r <- rank(u)
    stats::qnorm((r - 0.5) / length(r))
  }

  mat_norm <- apply(mat_std, 2, rank_norm)
  rownames(mat_norm) <- rownames(mat_std)

  # --- 6. Write Output ---
  # Match metadata to the kept features
  meta_final <- raw_subset[match(rownames(mat_norm), raw_subset$GeneID), required_cols, with = FALSE]

  # Format: GeneID, Chr, Start, End, Sample1, Sample2...
  out_df <- cbind(GeneID = meta_final$GeneID, as.data.frame(mat_norm))
  out_df_meta <- meta_final[, c("GeneID", "Chr", "Start", "End"), with = FALSE]

  data.table::fwrite(out_df, file = file_out, sep = "\t", quote = FALSE)
  data.table::fwrite(out_df_meta, file = file_out_meta,, sep = "\t", quote = FALSE)

  message("Success! Processed data written to: ", file_out)
  message("Success! Processed data meta written to: ", file_out_meta)
  invisible(out_df)
}
