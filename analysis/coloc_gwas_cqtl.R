##############################################
########### coloc between cQTL and GWAS hits of 36 traits ###########
##############################################
# load packages -----
# rm(list = ls())
library(tidyverse)


# I/O & paras -----
# ------------------------------------------------------------------
# Input files and parameters explanation
#
# --gwas_coloc_reg (file_gwas_coloc_reg):
#   GWAS regional summary statistics formatted for colocalization.
#   Each row corresponds to a SNP within a GWAS locus.
#
# --cwindow_fdr (file_cwindow_fdr):
#   Genome-wide cWindow (CACTI) QTL association results.
#   Each row represents a non-overlapping genomic window.
#
# --window_meta (file_window_meta):
#   RDS file containing metadata for cWindow groups.
#   Used to annotate genomic coordinates of each window.
#
# --p_qtl (file_p_qtl):
#   QTL SNP-level summary statistics.
#   Each row corresponds to a SNP tested in the QTL analysis.
#
# --af_qtl (file_af_qtl):
#   Allele frequency file for QTL SNPs.
#   Used to supply MAF information for colocalization.
#
# --chr:
#   Chromosome to analyze (e.g. "chr21").
#
# --gwasPhenocode:
#   GWAS trait identifier, including sample size information
#   (e.g. "mpv_N164454").
#
# --suffix:
#   String appended to output for bookkeeping (e.g. mark type).
#
# --qtlN:
#   Sample size of the QTL study.
#
# --fdr_level:
#   FDR threshold used to define significant cWindow QTLs.
#
# --dis_cis_qtl:
#   Cis window size for QTL–GWAS pairing.
#   Must be specified in kilobases (e.g. "500kb").
#
# --coloc_res (file_coloc_res):
#   Output file path for colocalization results.
#   Results include posterior probabilities (PP0–PP4)
#   for each GWAS–QTL locus pair.
# ------------------------------------------------------------------

if(interactive()){
  args <- scan(
    text = '
    ',
    what = 'character'
  )
} else{
  library("optparse")
  
  option_list <- list(
    make_option(c("--gwas_coloc_reg", "-g"), action = "store", type = "character", default = NULL,
                help = "gwas_coloc_reg"),
    make_option(c("--cwindow_fdr", "-w"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--window_meta"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--p_qtl"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--af_qtl"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--chr"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--gwasPhenocode"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--suffix"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--qtlN"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--fdr_level"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--dis_cis_qtl"), action = "store", type = "character", default = NULL,
                help = ""),
    make_option(c("--coloc_res"), action = "store", type = "character", default = NULL,
                help = "Output filename")
  )
  opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = FALSE)
}


file_gwas_coloc_reg <- opt$gwas_coloc_reg

file_cwindow_fdr <- opt$cwindow_fdr
file_window_meta <- opt$window_meta

file_p_qtl <- opt$p_qtl
file_af_qtl <- opt$af_qtl


chr <- str_extract(opt$chr, '\\d+') |> as.numeric()
gwasPhenocode <- opt$gwasPhenocode
suffix <- opt$suffix
qtlN <- as.numeric(opt$qtlN)
fdr_level <- as.numeric(opt$fdr_level)
dis_cis_qtl <- opt$dis_cis_qtl


## paras -----
dis_reg <- 1e+6
qtlType <- "quant"
gwasType <- "quant"


## output -----
file_coloc_res <- opt$coloc_res



# read files -----
# gwas reg files
gwas <- data.table::fread(file_gwas_coloc_reg, header = TRUE) %>%
  filter(chr == !!chr)

# qtl reg files
cwindow_fdr <- data.table::fread(file_cwindow_fdr)
window_meta <- readRDS(file_window_meta)

# cwindow
cwindow_fdr <- mutate(cwindow_fdr, if_sig = q <= fdr_level) %>%
  left_join(window_meta, by = c('group') )

# qtl sum stats
p_qtl <- data.table::fread(file_p_qtl, header = TRUE)

# af info
af_qtl <- data.table::fread(file_af_qtl, header = TRUE)

# add af info
p_qtl <- left_join(
  p_qtl, distinct(af_qtl, ID, CHROM, POS, AF),
  by = c('snp' = 'ID')
)

# add snp id
p_qtl$snp_simple <- paste(str_extract(p_qtl$CHROM, "\\d+"), p_qtl$POS, sep = ":")

# GWAS trait info, sample size
n_cases <- str_extract(gwasPhenocode, '_N\\d+') %>% str_extract('\\d+') %>% as.numeric()
n_controls <- NA
gwasN <- sum(n_cases, n_controls, na.rm = TRUE)

# re-format cis distance params
if(str_detect(dis_cis_qtl, '^\\d+[k, K][b, B]$')){
  dis_cis_qtl <- as.numeric(str_extract(dis_cis_qtl, '^\\d+')) * 1000
}else stop("Specify window size in kb. \n")




# find all pairs of loci and window that overlap -----
## left_join gwas regs, add qtl reg to gwas reg (1Mb)
## by distance of top gwas to top qtl
gwas_reg <- distinct(gwas, Region) %>%
  separate(
    col = Region, 
    into = c('reg_chr', 'reg_center'), 
    sep = ":", 
    remove = FALSE, convert = TRUE
  )

gwas_reg$qtl_reg <- mapply(
  FUN = function(x, y){
    filter(cwindow_fdr, phe_chr == str_glue("chr{x}")) %>%
      filter(
        abs(y - phe_from) < dis_cis_qtl + dis_reg/2 |
          abs(y - phe_to) < dis_cis_qtl + dis_reg/2
      ) %>%
      pull(group) %>%
      paste(collapse = ";")
  }, 
  gwas_reg$reg_chr, gwas_reg$reg_center
)

gwas_qtl_reg <- separate_rows(gwas_reg, qtl_reg, sep = ";")
n_region <- nrow(gwas_qtl_reg)



# coloc -----
## overlap gwas reg snps with qtl snps -----
gwas_qtl <- filter(gwas, (SNP_ID %in% !!p_qtl$snp_simple))
if(nrow(gwas_qtl) == 0) {
  cat(NULL, file_coloc_res)
  stop("No variants overlap between GWAS trait and QTL trait!")
}


# interate through regions of a chromosome
resColoc <- NULL
for(k_reg in 1:n_region){
  reg_gwas = gwas_qtl_reg$Region[k_reg]
  reg_qtl = gwas_qtl_reg$qtl_reg[k_reg]
  
  ## inner join regs for coloc (overlap snps between qtl and gwas) -----
  coloc_qtl_gwas = inner_join(
    filter(gwas_qtl, Region == !!reg_gwas),
    filter(p_qtl, group == !!reg_qtl),
    by = c('SNP_ID' = 'snp_simple')
  ) %>%
    distinct(SNP_ID, .keep_all = TRUE)
  if(nrow(coloc_qtl_gwas) == 0) next
  
  
  ## coloc -----
  D1 = list("pvalues" = replace(coloc_qtl_gwas$p_pco, coloc_qtl_gwas$p_pco == 0, min(coloc_qtl_gwas$p_pco[coloc_qtl_gwas$p_pco != 0])/10),
            "N" = qtlN,
            "MAF" = if_else(coloc_qtl_gwas$AF > 0.5, 1 - coloc_qtl_gwas$AF, coloc_qtl_gwas$AF),
            "type" = qtlType,
            "snp" = coloc_qtl_gwas$SNP_ID)
  D2 = list("pvalues" = coloc_qtl_gwas$pval,
            "beta" = coloc_qtl_gwas$beta,
            "varbeta" = (coloc_qtl_gwas$se)^2,
            "type" = gwasType,
            "s" = if(gwasType=="cc") n_cases/gwasN else NULL,
            "snp" = coloc_qtl_gwas$SNP_ID,
            "MAF" = coloc_qtl_gwas$maf,
            "N" = gwasN)
  
  if(gwasType=="cc"){
    D2 = within(D2, rm(MAF))
  }else{
    D2 = within(D2, rm(s))
  }
  
  # do coloc
  coloc_res = coloc::coloc.abf(D1, D2)
  
  # progress
  cat(str_glue('On chr{chr}, {k_reg}-th reg (out of {n_region} regions) is done... \n\n'))
  
  
  ## aggregate results across regions -----
  resColoc = c(
    resColoc,
    list(
      enframe(c(
        'reg_qtl' = reg_qtl, 'Region' = reg_gwas, 
        'chr' = chr, 
        'minp_qtl_coloc' = min(coloc_qtl_gwas$p_pco), 'minp_gwas_coloc' = min(coloc_qtl_gwas$pval),
        coloc_res$summary
      ))
    )
  )
}



# save -----
bind_rows(resColoc) %>%
  group_by(name) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = name, values_from = value) %>%
  select(-id) %>%
  mutate('Phenocode' = gwasPhenocode, 'suffix' = suffix) %>%
  data.table::fwrite(file = file_coloc_res, quote = FALSE, sep = "\t")

