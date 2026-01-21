#' Run Nominal Cis-QTL Mapping (CACTI-S)
#'
#' Performs linear regression using MatrixEQTL to generate summary statistics.
#'
#' @param file_pheno Path to the processed phenotype file with position meta and expression data (from `cacti_s_preprocess`).
#' @param file_pheno_meta Path for the meta file of processed phenotype.
#' @param file_cov Path to covariate matrix.
#' @param file_vcf Path to input VCF file.
#' @param file_qtl_out Output path for summary statistics.
#' @param file_geno (Optional) Path to genotype matrix if no VCF.
#' @param file_snp_pos (Optional) Path to SNP positions if no VCF.
#' @param cis_dist Cis-window distance (default 100000 bp = 100kb).
#' @param p_threshold P-value threshold for output (default 1.0, print all associations).
#'
#' @return Invisibly returns the summary stats.
#' @export
#'
#' @examples
#' \dontrun{
#' # Use local paths for testing
#' file_pheno <- "inst/extdata/test_results/test_pheno_norm.txt"
#' file_pheno_meta <- "inst/extdata/test_results/test_pheno_norm_meta.txt"
#' file_cov   <- "inst/extdata/test_cov.txt"
#' file_vcf   <- "inst/extdata/test_geno.vcf"
#'
#' file_out <- "inst/extdata/test_results/test_cis_qtl_stats.txt"
#'
#' # Run
#' if (file.exists(file_pheno) && file.exists(file_vcf)) {
#'   cacti_s_map_cis(
#'     file_pheno   = file_pheno,
#'     file_pheno_meta = file_pheno_meta,
#'     file_cov     = file_cov,
#'     file_vcf     = file_vcf,
#'     file_qtl_out = file_out,
#'     cis_dist     = 100000
#'   )
#' }
#' }
cacti_s_map_cis <- function(
    file_pheno,
    file_pheno_meta,
    file_cov,
    file_qtl_out,
    file_vcf = NULL,
    file_geno = NULL,
    file_snp_pos = NULL,
    cis_dist = 100000,
    p_threshold = 1.0
) {
  if (!requireNamespace("MatrixEQTL", quietly = TRUE)) stop("MatrixEQTL required.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required.")

  # --- 1. Prepare Phenotype Data ---
  message("Preparing phenotype data...")

  # Skip meta cols, assume standard output from cacti_s_preprocess: GeneID (1) and columns 5+
  dt_pheno <- data.table::fread(file_pheno, header = TRUE)
  file_pheno_clean <- tempfile(fileext = ".txt")
  data.table::fwrite(dt_pheno, file_pheno_clean, sep = "\t", quote = FALSE)

  target_samples <- colnames(dt_pheno)[2:ncol(dt_pheno)]
  message("  Phenotype samples: ", length(target_samples))

  # --- 2. Handle Genotypes (VCF or Text) ---
  clean_geno <- FALSE; clean_pos <- FALSE
  if (!is.null(file_vcf)) {
    if (file_vcf == "" || !file.exists(file_vcf)) stop("VCF file not found: ", file_vcf)
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) stop("VariantAnnotation required.")
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required.")
    if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) stop("GenomeInfoDb required.")
    if (!requireNamespace("BiocGenerics", quietly = TRUE)) stop("BiocGenerics required.")

    message("Converting VCF to genotype matrix...")
    file_geno <- tempfile(fileext = ".txt")
    file_snp_pos <- tempfile(fileext = ".txt")
    clean_geno <- TRUE; clean_pos <- TRUE

    # Read VCF
    vcf <- VariantAnnotation::readVcf(file_vcf)

    # Convert to numeric (0,1,2)
    gts_list <- VariantAnnotation::genotypeToSnpMatrix(vcf)
    mat_numeric <- as(gts_list$genotypes, "numeric")

    # Sample matching
    vcf_samples <- rownames(mat_numeric)
    common_samples <- intersect(target_samples, vcf_samples)

    if (length(common_samples) == 0) stop("No overlapping samples between Phenotype and VCF.")

    mat_numeric <- mat_numeric[common_samples, , drop = FALSE]
    mat_geno <- t(mat_numeric) # SNPs x Samples

    # Ensure SNP positions match the SNPs in mat_geno
    rr <- SummarizedExperiment::rowRanges(vcf)
    keep_snps <- rownames(mat_geno)
    rr_sub <- rr[keep_snps, ]

    if (length(keep_snps) != length(rr_sub)) {
      warning("Mismatch between genotype matrix rows and VCF rowRanges. Trying to align by name...")
      rr_sub <- rr[match(keep_snps, names(rr)), ]
    }

    # Create Genotype File
    df_geno <- data.frame(ID = rownames(mat_geno), mat_geno, check.names = FALSE)
    data.table::fwrite(df_geno, file_geno, sep = "\t", quote = FALSE)

    # Create SNP Position File
    df_pos <- data.frame(snp = names(rr_sub),chr = as.character(GenomeInfoDb::seqnames(rr_sub)),pos = BiocGenerics::start(rr_sub))
    data.table::fwrite(df_pos, file_snp_pos, sep = "\t", quote = FALSE)

    # Subset phenotype file to common samples with VCF
    if (length(common_samples) < length(target_samples)) {
      warning("Subsetting phenotype to ", length(common_samples), " common samples.")
      dt_pheno_sub <- dt_pheno[, c("GeneID", common_samples), with = FALSE]
      data.table::fwrite(dt_pheno_sub, file_pheno_clean, sep = "\t", quote = FALSE)
    }
  } else {
    if (is.null(file_geno)) stop("Must provide file_vcf or file_geno.")
  }


  # --- 3. Handle Covariates ---
  # MatrixEQTL crashes if Covariates have different samples than Pheno/Geno.
  file_cov_clean <- NULL
  if (!is.null(file_cov) && file.exists(file_cov)) {
    message("Preparing covariates...")
    cov_dt <- data.table::fread(file_cov)

    # Column 1 is ID
    cov_id_col <- names(cov_dt)[1]
    cov_samples <- names(cov_dt)[-1]

    # Check coverage
    missing_cov <- setdiff(common_samples, cov_samples)
    if (length(missing_cov) > 0) {
      stop("Covariate file is missing samples present in Pheno/VCF: ", paste(head(missing_cov), collapse=", "))
    }

    # Subset columns: ID + Common Samples (in exact order)
    cov_subset <- cov_dt[, c(cov_id_col, common_samples), with = FALSE]

    file_cov_clean <- tempfile(fileext = ".txt")
    data.table::fwrite(cov_subset, file_cov_clean, sep = "\t", quote = FALSE)
  }


  # --- 4. Load Data into MatrixEQTL ---

  # A. Genotypes
  # Format: ID, S1, S2...
  snps <- MatrixEQTL::SlicedData$new()
  snps$fileDelimiter <- "\t"
  snps$fileOmitCharacters <- "NA"
  snps$fileSkipRows <- 1
  snps$fileSkipColumns <- 1
  snps$fileSliceSize <- 2000
  snps$LoadFile(file_geno)

  # B. Phenotypes
  # Format: ID, S1, S2...
  gene <- MatrixEQTL::SlicedData$new()
  gene$fileDelimiter <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows <- 1
  gene$fileSkipColumns <- 1
  gene$fileSliceSize <- 2000
  gene$LoadFile(file_pheno_clean)

  # C. Covariates
  # Format: ID, S1, S2...
  cvrt <- MatrixEQTL::SlicedData$new()
  cvrt$fileDelimiter <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows <- 1
  cvrt$fileSkipColumns <- 1
  if (!is.null(file_cov_clean)) cvrt$LoadFile(file_cov_clean)


  # D. Positions
  snpspos <- data.table::fread(file_snp_pos, header = TRUE, data.table = FALSE)
  colnames(snpspos) <- c("snp", "chr", "pos")

  # Extract gene positions from ORIGINAL pheno file (Cols 1-4)
  genepos <- data.table::fread(file_pheno_meta, header = TRUE, data.table = FALSE)
  colnames(genepos) <- c("geneid", "chr", "left", "right")

  # --- 4. Run MatrixEQTL ---
  tmp_out <- tempfile()
  message("Running MatrixEQTL...")

  me <- MatrixEQTL::Matrix_eQTL_main(
    snps = snps, gene = gene, cvrt = cvrt,
    output_file_name = NULL, pvOutputThreshold = 0,
    useModel = MatrixEQTL::modelLINEAR, errorCovariance = numeric(), verbose = TRUE,
    output_file_name.cis = tmp_out, pvOutputThreshold.cis = p_threshold,
    snpspos = snpspos, genepos = genepos, cisDist = cis_dist,
    pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE
  )

  # --- 5. Format Output ---
  res <- data.table::fread(tmp_out)
  out_df <- res |> dplyr::select(phe_id = gene, var_id = SNP, z = `t-stat`, pval = `p-value`)
  data.table::fwrite(out_df, file = file_qtl_out, sep = "\t", quote = FALSE)

  # Clean up
  unlink(tmp_out); unlink(file_pheno_clean)
  if (clean_geno) unlink(file_geno)
  if (clean_pos) unlink(file_snp_pos)

  message("Success! Stats written to: ", file_qtl_out)

  invisible(out_df)
}
