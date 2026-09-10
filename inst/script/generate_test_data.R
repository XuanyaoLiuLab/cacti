# ==============================================================================
# Script to Generate Bundled Test Data for CACTI-S
# ==============================================================================
# Dependencies: Rsamtools
# ==============================================================================

# Setup Output Directory
# -------------------------
out_dir <- "inst/extdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Generating test data in: ", out_dir)


# -------------------------------------------------------------
# 1. Generate BAMs
# -------------------------------------------------------------
if (requireNamespace("Rsamtools", quietly = TRUE)) {
  # Twelve genuine 5-kb segments, with a different profile for each sample.
  # The former duplicated sample pairs and two retained segments produced
  # exact genotype-phenotype fits after normalization (infinite t-statistics).
  # This synthetic fixture checks execution, not statistical calibration.
  saf_df <- data.frame(
    Start = seq(1L, by = 5000L, length.out = 12L),
    End = seq(5000L, by = 5000L, length.out = 12L)
  )
  set.seed(20260909)
  toy_counts <- matrix(sample(10:90, 12L * 4L, replace = TRUE), nrow = 12L)

  # -------------------------
  create_fake_bam <- function(filename, saf, n_reads_vec) {
    header <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:100000")
    alignments <- c()
    for (i in 1:nrow(saf)) {
      if (n_reads_vec[i] > 0) {
        for (r in 1:n_reads_vec[i]) {
          line <- paste(paste0("seg", i, "_read", r), 0, "chr1",
                        saf$Start[i] + 10L + r, 60, "50M", "*", 0, 0,
                        paste(rep("A", 50), collapse = ""),
                        paste(rep("I", 50), collapse = ""), sep = "\t")
          alignments <- c(alignments, line)
        }
      }
    }
    sam_file <- paste0(filename, ".sam")
    bam_dest <- gsub(".bam$", "", filename)

    if (file.exists(sam_file)) unlink(sam_file)
    if (file.exists(paste0(bam_dest, ".bam"))) unlink(paste0(bam_dest, ".bam"))
    if (file.exists(paste0(bam_dest, ".bam.bai"))) unlink(paste0(bam_dest, ".bam.bai"))

    # Write SAM
    writeLines(c(header, alignments), sam_file)

    # Convert to BAM
    Rsamtools::asBam(sam_file, destination = bam_dest, overwrite = TRUE, indexDestination = TRUE)
    unlink(sam_file)
  }

  # Generate 4 Samples
  # ----------------------
  for (j in seq_len(ncol(toy_counts))) {
    create_fake_bam(file.path(out_dir, paste0("Sample", j, ".bam")),
                    saf_df, toy_counts[, j])
  }

  message("[Part A] Four BAMs and indexes created across twelve 5-kb segments.")
} else {
  warning("Rsamtools not installed. Skipping BAM generation.")
}



# -------------------------------------------------------------
# 2. Generate VCF (Genotypes)
# -------------------------------------------------------------
# Setup 50 Samples
n_samples <- 50
samples <- paste0("Sample", 1:n_samples)


# 50 SNPs on chr1
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste(samples, collapse = "\t"))
)

# Simulate genotypes (0/0, 0/1, 1/1)
# 50 SNPs at pos 10,000 to 60,000
pos <- seq(10000, 60000, by = 1000)
vcf_lines <- character(length(pos))

set.seed(42)
for (i in seq_along(pos)) {
  gts <- sample(c("0/0", "0/1", "1/1"), n_samples, replace = TRUE)
  vcf_lines[i] <- sprintf(
    "chr1\t%d\tsnp_%d\tA\tG\t.\tPASS\t.\tGT\t%s",
    pos[i], i, paste(gts, collapse = "\t")
  )
}
writeLines(c(vcf_header, vcf_lines), file.path(out_dir, "test_geno.vcf"))


# -------------------------------------------------------------
# 3. Generate Covariates
# -------------------------------------------------------------
# Covariates (PC1)
cov_df <- data.frame(ID = "PC1", matrix(rnorm(1 * n_samples), nrow = 1))
colnames(cov_df) <- c("ID", samples)
write.table(cov_df, file.path(out_dir, "test_cov.txt"), sep="\t", quote=FALSE, row.names=FALSE)

message("[Part B] Created VCF, Pheno, and Covariates.")
message("All test data generated in 'inst/extdata'.")


# -------------------------------------------------------------
# 4. Generate Fake Peak BED File For Filtering Mode `match_overlap_count`
# -------------------------------------------------------------
n_peaks = 5
chr = "chr1"

starts <- seq(1000, by = 2000, length.out = n_peaks)
ends   <- starts + 500

df_bed <- data.frame(chr   = rep(chr, n_peaks), start = starts, end   = ends)

write.table(df_bed, file = file.path(out_dir, "test_peaks.bed"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
