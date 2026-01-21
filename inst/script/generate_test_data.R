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
  # -------------------------
  create_fake_bam <- function(filename, saf, n_reads_vec) {
    header <- c("@HD\tVN:1.0\tSO:coordinate", "@SQ\tSN:chr1\tLN:100000")
    alignments <- c()
    for (i in 1:nrow(saf)) {
      if (n_reads_vec[i] > 0) {
        for (r in 1:n_reads_vec[i]) {
          line <- sprintf("read_%d_%d\t0\tchr1\t%d\t255\t50M\t*\t0\t0\t*\t*", i, r, saf$Start[i] + 10)
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
    Rsamtools::asBam(sam_file, destination = bam_dest, overwrite = TRUE)
    unlink(sam_file)
  }

  # Generate 4 Samples
  # ----------------------
  # --- FIX: Rename BAMs to match VCF sample IDs (Sample1, Sample2) ---
  create_fake_bam(file.path(out_dir, "Sample1.bam"), saf_df, c(10, 20, 5, 0, 15))
  create_fake_bam(file.path(out_dir, "Sample2.bam"), saf_df, c(5, 50, 2, 10, 8))
  create_fake_bam(file.path(out_dir, "Sample3.bam"), saf_df, c(10, 20, 5, 0, 15))
  create_fake_bam(file.path(out_dir, "Sample4.bam"), saf_df, c(5, 50, 2, 10, 8))

  message("[Part A] BAMs and SAF created.")
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
