#' Run MatrixEQTL and write CACTI-compatible cis-QTL summary stats
#'
#' Thin wrapper around [cacti_s_map_cis()] that generates a summary-statistics
#' table with columns `phe_id`, `var_id`, `z`, `pval` for downstream CACTI
#' peak-window testing.
#'
#' @inheritParams cacti_s_map_cis
#'
#' @return Invisibly returns the summary-statistics data frame.
#' @export
cacti_matrixqtl_cis <- function(
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
  cacti_s_map_cis(
    file_pheno = file_pheno,
    file_pheno_meta = file_pheno_meta,
    file_cov = file_cov,
    file_qtl_out = file_qtl_out,
    file_vcf = file_vcf,
    file_geno = file_geno,
    file_snp_pos = file_snp_pos,
    cis_dist = cis_dist,
    p_threshold = p_threshold
  )
}
