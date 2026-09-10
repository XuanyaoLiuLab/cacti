# Regenerate bundled example results and the complete documentation website.
# Run from the package source directory: Rscript --vanilla build_pkg.R
# Requires the package dependencies plus pkgload, pkgdown, knitr, rmarkdown and Pandoc.
# This does not install dependencies, delete old files, run CI, or push to GitHub.

stopifnot(file.exists("DESCRIPTION"), read.dcf("DESCRIPTION")[1, "Package"] == "cacti")
for (package in c("pkgload", "pkgdown", "knitr", "rmarkdown")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install build dependency: ", package)
}
if (!rmarkdown::pandoc_available()) stop("Pandoc is required to build the website.")
pkgload::load_all(".", export_all = FALSE)

fixture <- function(name) {
  path <- file.path("inst", "extdata", name)
  if (!file.exists(path)) stop("Missing bundled input: ", path)
  normalizePath(path, mustWork = TRUE)
}
read_output <- function(path, columns) {
  stopifnot(length(path) == 1L, file.exists(path), file.info(path)$size > 0)
  result <- data.table::fread(path, data.table = FALSE)
  stopifnot(nrow(result) > 0L, all(columns %in% names(result)))
  result
}
check_probability <- function(x) {
  stopifnot(is.numeric(x), length(x) > 0L, all(is.finite(x)), all(x >= 0 & x <= 1))
}

out_dir <- file.path("inst", "extdata", "test_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. CACTI-S: regenerate counts, normalized phenotypes, and single-segment statistics.
segment_prep <- cacti::cacti_s_prep(
  file_bams = vapply(paste0("Sample", 1:4, ".bam"), fixture, character(1)),
  file_vcf = fixture("test_geno.vcf"),
  file_cov = fixture("test_cov.txt"),
  out_dir = out_dir,
  out_prefix = "test",
  genome = "hg19",
  segment_size = "5kb",
  threads = 1
)
raw_counts <- read_output(segment_prep$raw_counts,
                         c("GeneID", "Chr", "Start", "End", "Strand", "Length"))
count_columns <- setdiff(names(raw_counts), c("GeneID", "Chr", "Start", "End", "Strand", "Length"))
count_matrix <- as.matrix(raw_counts[count_columns])
stopifnot(ncol(count_matrix) == 4L, all(is.finite(count_matrix)), all(count_matrix >= 0))
nonzero_segments <- raw_counts$GeneID[rowSums(count_matrix) > 0]
phenotype <- read_output(segment_prep$pheno, c("GeneID", paste0("Sample", 1:4)))
stopifnot(length(nonzero_segments) == 12L,
          setequal(phenotype$GeneID, nonzero_segments),
          all(is.finite(as.matrix(phenotype[-1]))))
segment_stats <- read_output(segment_prep$qtl_stats, c("phe_id", "var_id", "z", "pval"))
check_probability(segment_stats$pval)
stopifnot(all(is.finite(segment_stats$z)),
          max(abs(2 * pnorm(-abs(segment_stats$z)) - segment_stats$pval)) < 1e-12,
          setequal(segment_stats$phe_id, phenotype$GeneID),
          !anyDuplicated(segment_stats[c("phe_id", "var_id")]))

# 2. Group adjacent 5-kb segments into 10-kb windows for a genuine multisegment test.
segment_windows <- cacti::cacti_run_chr(
  window_size = "10kb",
  file_pheno_meta = file.path(out_dir, "test_pheno_norm_meta.txt"),
  file_pheno = segment_prep$pheno,
  file_cov = fixture("test_cov.txt"),
  qtl_file = segment_prep$qtl_stats,
  chr = "chr1",
  min_peaks = 2,
  out_prefix = file.path(out_dir, "cactis_chr1")
)
groups <- read_output(segment_windows$file_peak_group, c("group", "n_pid_group"))
window_stats <- read_output(segment_windows$file_p_peak_group, c("group", "snp", "pval"))
check_probability(window_stats$pval)
stopifnot(all(groups$n_pid_group == 2L),
          all(window_stats$pval > 0),
          sum(groups$n_pid_group) == nrow(phenotype),
          setequal(window_stats$group, groups$group),
          !anyDuplicated(window_stats[c("group", "snp")]))
message("CACTI-S regenerated: ", nrow(phenotype), " segments, ", nrow(segment_stats),
        " SNP-segment associations, ", nrow(window_stats), " SNP-window associations in ",
        nrow(groups), " windows.")

# 3. CACTI: regenerate the bundled signed-Z input and its window results with FDR.
cacti::cacti_matrixqtl_cis(
  file_pheno_meta = fixture("test_cacti_peak_chr5_pheno_meta.bed"),
  file_pheno = fixture("test_cacti_peak_chr5_pheno.txt"),
  file_cov = fixture("test_cacti_peak_chr5_covariates.txt"),
  file_vcf = fixture("test_cacti_peak_chr5_geno.vcf"),
  file_qtl_out = file.path("inst", "extdata", "test_cacti_peak_chr5_matrixqtl_sumstats.txt.gz"),
  cis_dist = 100000,
  p_threshold = 1.0
)
peak_results <- cacti::cacti_peak_window(
  window_size = "50kb",
  file_pheno_meta = fixture("test_cacti_peak_chr5_pheno_meta.bed"),
  file_pheno = fixture("test_cacti_peak_chr5_pheno.txt"),
  file_cov = fixture("test_cacti_peak_chr5_covariates.txt"),
  qtl_file = fixture("test_cacti_peak_chr5_matrixqtl_sumstats.txt.gz"),
  chr = "chr5",
  do_fdr = TRUE,
  out_prefix = file.path(out_dir, "test_cacti_peak_chr5_from_sumstats")
)
peak_stats <- read_output(peak_results$file_p_peak_group, c("group", "snp", "pval"))
peak_fdr <- read_output(peak_results$file_fdr_out, c("group", "best_hit_acat", "q"))
check_probability(peak_stats$pval)
check_probability(peak_fdr$best_hit_acat)
check_probability(peak_fdr$q)
stopifnot(setequal(peak_fdr$group, peak_stats$group))
message("CACTI regenerated: ", nrow(peak_stats), " SNP-window associations in ",
        nrow(peak_fdr), " windows, with FDR.")

# 4. Rebuild both articles and the rest of the site through the same pkgdown template.
# Install the current source into pkgdown's temporary library, and evaluate the
# articles in fresh processes. Do not substitute saved output or reuse cached pages.
pkgdown::build_site(lazy = FALSE, quiet = FALSE, preview = FALSE,
                    new_process = TRUE, install = TRUE)
