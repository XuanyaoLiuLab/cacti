#' Run CACTI-S Prep Pipeline
#'
#' A master wrapper function that executes the full CACTI-S preparation pipeline.
#' It automates the transition from raw BAM files to a normalized phenotype matrix
#' and optionally performs nominal cis-QTL mapping.
#'
#' **Pipeline Stages:**
#' \enumerate{
#'   \item \strong{Create Segments:} Generates fixed-size genomic windows (SAF) for the specified genome build.
#'   \item \strong{Count Reads:} Quantifies reads in these segments using \code{Rsubread::featureCounts}.
#'   \item \strong{Preprocess:} Performs QC (removing zero-count segments), Normalization (CPM),
#'     Filtering (different modes), Standardization (Z-score), and Rank Normalization (INT).
#'   \item \strong{Map Cis-QTLs:} (Optional) Runs linear regression using \code{MatrixEQTL} if VCF and Covariates are provided.
#' }
#'
#' @param file_bams Character vector of file paths to the input BAM files. These should be sorted and indexed.
#' @param out_dir Directory path where all output files will be saved. Created automatically if it does not exist.
#' @param out_prefix String prefix for all generated filenames (default "cacti_s").
#'   Files created: \code{_segments.saf}, \code{_raw_counts.txt}, \code{_pheno_norm.txt}, \code{_cis_qtl_stats.txt}, etc.
#' @param file_saf_custom (Optional) Path to an existing SAF file. If provided, segment generation is skipped.
#' @param file_vcf (Optional) Path to the input VCF file containing genotype data.
#'   Required if you want to run the mapping step.
#' @param file_cov (Optional) Path to the covariate matrix file.
#'   Format: Rows are covariates (ID, S1, S2...), Columns are samples. Required if \code{file_vcf} is provided.
#' @param genome Genome build to use for segment generation. Options: \code{"hg19"} (default) or \code{"hg38"}.
#' @param segment_size String defining the window size for tiling the genome (e.g., \code{"1kb"}, \code{"5kb"}, \code{"10kb"}).
#'   Default is \code{"5kb"}.
#' @param filter_mode Strategy for selecting informative segments during preprocessing:
#'   \itemize{
#'     \item \code{"cpm_threshold"} (Default): Keeps segments with Counts Per Million (CPM) > \code{min_cpm}
#'       in at least \code{min_prop} proportion of samples.
#'     \item \code{"match_overlap_count"}: Calculates $N$, the number of segments that physically overlap
#'       with regions in \code{filter_bed}. It then selects the top $N$ segments based on average CPM intensity.
#'     \item \code{"prop"}: Simply keeps the top \code{prop_top} proportion of segments ranked by average CPM.
#'   }
#' @param filter_bed (Required only for \code{match_overlap_count}) Path to a BED file (e.g., ATAC-seq peaks,
#'   H3K27ac peaks, or known QTL regions) used to determine the number of segments to keep.
#' @param min_cpm (For \code{cpm_threshold}) Minimum CPM intensity threshold. Default 1.
#' @param min_prop (For \code{cpm_threshold}) Minimum proportion of samples (0-1) that must exceed \code{min_cpm}. Default 0.2.
#' @param prop_top (For \code{prop}) Proportion (0-1) of top-expressed segments to keep. Default 0.2.
#' @param cis_dist (For Mapping) The maximum distance (bp) between a gene (segment) and a SNP to be considered "cis". Default 100,000 bp.
#' @param threads Integer, number of threads to use for the read counting step. Default 1.
#' @param isPairedEnd Logical, indicating if the BAM files contain paired-end reads. Default \code{FALSE}.
#'
#' @return Invisibly returns a list containing paths to the generated files:
#'   \itemize{
#'     \item \code{segments.saf}: Path to the segment definition file.
#'     \item \code{raw_counts.txt}: Path to the featureCounts output.
#'     \item \code{pheno_norm.txt}: Path to the final normalized phenotype matrix.
#'     \item \code{cis_qtl_stats.txt}: (If mapped) Path to the cis-QTL summary statistics.
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Locate test data bundled with the package
#' file_bam1 <- system.file("extdata", "Sample1.bam", package = "cacti")
#' file_bam2 <- system.file("extdata", "Sample2.bam", package = "cacti")
#' file_bam3 <- system.file("extdata", "Sample3.bam", package = "cacti")
#' file_bam4 <- system.file("extdata", "Sample4.bam", package = "cacti")
#' file_vcf  <- system.file("extdata", "test_geno.vcf", package = "cacti")
#' file_cov  <- system.file("extdata", "test_cov.txt", package = "cacti")
#'
#' # Create a temporary output directory
#' out_dir <- "inst/extdata/test_results"
#'
#' if (nchar(file_bam1) > 0) {
#'   # --- Run Full Pipeline (BAM -> QTL Stats) ---
#'   res <- cacti_s_prep(
#'     file_bams = c(file_bam1, file_bam2, file_bam3, file_bam4),
#'     file_vcf = file_vcf,
#'     file_cov = file_cov,
#'
#'     out_dir = out_dir,
#'     out_prefix = "test",
#'
#'     # Segment Parameters
#'     file_saf_custom = NULL,
#'     genome = "hg19",
#'     segment_size = "5kb",
#'
#'     # Filtering Parameters
#'     filter_mode = "cpm_threshold",
#'     min_cpm = 1,
#'     min_prop = 0.2,
#'
#'     # Execution Parameters
#'     threads = 1
#'   )
#' }
#' }
cacti_s_prep <- function(
    file_bams,
    file_vcf = NULL,
    file_cov = NULL,
    out_dir,
    out_prefix = "cacti_s",
    file_saf_custom = NULL,
    genome = "hg19",
    segment_size = "5kb",
    filter_mode = "cpm_threshold",
    filter_bed = NULL,
    min_cpm = 1,
    min_prop = 0.2,
    prop_top = 0.2,
    cis_dist = 100000,
    threads = 1,
    isPairedEnd = FALSE
) {
  # --- 0. Setup ---
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("Created output directory: ", out_dir)
  }

  # Define filenames
  file_saf     <- file.path(out_dir, paste0(out_prefix, "_segments.saf"))
  file_counts  <- file.path(out_dir, paste0(out_prefix, "_raw_counts.txt"))
  file_pheno   <- file.path(out_dir, paste0(out_prefix, "_pheno_norm.txt"))
  file_pheno_meta   <- file.path(out_dir, paste0(out_prefix, "_pheno_norm_meta.txt"))
  file_qtl     <- file.path(out_dir, paste0(out_prefix, "_cis_qtl_stats.txt"))

  message("=======================================================")
  message(" CACTI-S Preparation Pipeline")
  message("=======================================================")
  message("Output Directory: ", out_dir)
  message("Input BAMs: ", length(file_bams))

  # --- 1. Create Segments (or use custom) ---
  if (!is.null(file_saf_custom)) {
    message("\n[Step 1] Using custom segments: ", file_saf_custom)
    file_saf <- file_saf_custom
  } else {
    message("\n[Step 1] Creating Genomic Segments...")
    cacti_s_create_segments(
      file_saf_out = file_saf,
      genome = genome,
      segment_size = segment_size
    )
  }

  # --- 2. Count Reads ---
  message("\n[Step 2] Counting Reads...")
  cacti_s_count(
    file_bams = file_bams,
    file_saf = file_saf,
    file_counts_out = file_counts,
    threads = threads,
    isPairedEnd = isPairedEnd
  )

  # --- 3. Preprocess (QC/Norm/Filter) ---
  message("\n[Step 3] Preprocessing (QC, Filter, Normalize)...")
  cacti_s_preprocess(
    file_raw_counts = file_counts,
    file_out = file_pheno,
    file_out_meta = file_pheno_meta,
    filter_mode = filter_mode,
    filter_bed = filter_bed,
    min_cpm = min_cpm,
    min_prop = min_prop,
    prop_top = prop_top
  )

  # --- 4. Map Cis-QTLs (Optional) ---
  run_mapping <- !is.null(file_vcf) && !is.null(file_cov)
  if (run_mapping) {
    message("\n[Step 4] Mapping Cis-QTLs...")
    cacti_s_map_cis(
      file_pheno = file_pheno,
      file_pheno_meta = file_pheno_meta,
      file_cov = file_cov,
      file_vcf = file_vcf,
      file_qtl_out = file_qtl,
      cis_dist = cis_dist
    )
    message("Mapping complete: ", file_qtl)
  } else {
    message("\n[Step 4] Skipping Mapping (VCF or Covariates missing).")
  }


  # --- FINAL SUMMARY SECTION ---
  message("\n=======================================================")
  message(" Pipeline Complete!")
  message("=======================================================")
  message("Summary of Generated Files:")
  message("  [1] Segments (SAF):      ", file_saf)
  message("  [2] Raw Counts:          ", file_counts)
  message("  [3] Normalized Pheno:    ", file_pheno)
  if (run_mapping) {
    message("  [4] Cis-QTL Stats:       ", file_qtl)
  }
  message("=======================================================")


  # Prepare Output List
  out_list <- list(
    saf = file_saf,
    raw_counts = file_counts,
    pheno = file_pheno
  )
  if (run_mapping) {
    out_list$qtl_stats <- file_qtl
  }

  invisible(out_list)
}
