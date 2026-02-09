##############################################
########### Define gwas region for coloc ###########
##############################################
# load packages -----
# rm(list = ls())
library(tidyverse)

# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    baso_N171846_narrow_form.tsv.gz
    gwas_baso_N171846_snp_converted_hg38.bed
    reg_gwas_baso_N171846.txt.gz
    reg_gwas_baso_N171846_converted_hg38.txt.gz
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_gwas <- args[1]
file_conversion <- args[2]

## paras -----
p_included_thre <- 1e-5
dis_reg <- 1e+6

## output -----
file_gwas_coloc_reg <- args[3]
file_gwas_coloc_reg_liftover <- args[4]



# read files -----
# liftover file
conversion <- data.table::fread(
  file_conversion, 
  header = FALSE, 
  col.names = c('chr', 'pos_l', 'pos_converted', 'snp')
)

# gwas files
gwas_col <- c("VARIANT", "CHR", "BP", "MA_FREQ", "EFFECT", "SE", "P", "ID_dbSNP49")
gwas <- data.table::fread(
  cmd = paste("gunzip -c", file_gwas), 
  select = gwas_col
)
colnames(gwas) <- c("snp", "chr", "pos", "maf", "beta", "se", "pval", "rsid")
gwas$pval <- as.numeric(gwas$pval)



# 1. Prep GWAS trait -----
gwas$SNP_ID <- paste(gwas$chr, gwas$pos, sep = ":")

## remove bad snps -----
## gwas good SNPs & remove duplicated SNPs & remove SNPs whose se(\beta) equals 0, as coloc needs to use 1/se(\beta)
gwas <- filter(
  gwas,
  !is.na(pval) & !duplicated(SNP_ID) & !(se == 0)
)

## remove MHC region -----
## (GRCh37, chr6:28,477,797-33,448,354)
## (GRCh38.p14, chr6:28,510,120-33,480,577)
gwas <- filter(gwas, !(chr == 6 & between(pos, 28477797, 33448354)))

if(nrow(gwas) == 0) stop("No variants left for this GWAS trait!")



# 2. Define GWAS regions -----
## Find gwas lead SNPs and define coloc regions -----
tmpgwas <- arrange(gwas, pval)
gwas_coloc_reg <- NULL
k <- 0
lSignal <- tmpgwas[1, ]
while (lSignal$pval <= p_included_thre) {
  indReg <- tmpgwas$chr == lSignal$chr & 
    abs(tmpgwas$pos - lSignal$pos) < dis_reg/2
  
  gwas_coloc_reg <- c(
    gwas_coloc_reg, 
    list(
      filter(tmpgwas, !!indReg) %>% mutate(Region = !!lSignal$SNP_ID)
    )
  )
  
  tmpgwas <- tmpgwas[!indReg, ]
  
  lSignal <- tmpgwas[1, ]
  
  k <- k + 1
  cat(str_glue('{k}-th region done. P is {lSignal$pval}. \n\n'))
}
gwas_coloc_reg <- bind_rows(gwas_coloc_reg)


# 3. liftover regions -----
ind_reg_to_snp <- match(gwas_coloc_reg$Region, gwas_coloc_reg$SNP_ID)

left_join(
  gwas_coloc_reg,
  distinct(conversion, snp, pos_converted), 
  by = c('snp')
) %>%
  rename(snp_old = snp) %>%
  mutate(pos = pos_converted, pos_converted = NULL) %>%
  mutate(SNP_ID = paste(chr, pos, sep = ":")) %>%
  mutate(Region = SNP_ID[!!ind_reg_to_snp]) %>%
  filter(!is.na(pos)) %>%
  data.table::fwrite(file_gwas_coloc_reg_liftover, quote = FALSE, sep = "\t")


# 4. print out key message or write out -----
data.table::fwrite(gwas_coloc_reg, file_gwas_coloc_reg, quote = FALSE, sep = "\t")

