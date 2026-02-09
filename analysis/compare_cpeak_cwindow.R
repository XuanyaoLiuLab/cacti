##############################################
########### Compare cwindow and cpeak overlaps ###########
##############################################
# rm(list = ls())
# load packages -----
library(tidyverse)

# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_cwindow_fdr <- args[1] # cWindow association results. Each row corresponds to a non-overlapping genomic window.
file_cpeak_fdr <- args[2] # cPeak association results. Each row corresponds to a chromatin peak.
file_window_meta <- args[3] # File containing metadata for cWindow groups.
fdr_level <- as.numeric(args[4]) # FDR threshold used to define significant cWindows.

# read files -----
cwindow_fdr <- data.table::fread(file_cwindow_fdr)
cpeak_fdr <- data.table::fread(file_cpeak_fdr)
window_meta <- readRDS(file_window_meta)

# add cols -----
## 1. add cols of if_sig cpeaks / cwindows and phe center -----
## cwindow
cwindow_fdr <- mutate(cwindow_fdr, if_sig = q <= fdr_level) %>%
  left_join(window_meta, by = c('group') ) %>%
  rowwise() %>%
  mutate(phe_center = round(mean(c(phe_from, phe_to)))) %>%
  ungroup()

## cpeak
cpeak_fdr <- mutate(cpeak_fdr, if_sig = q <= fdr_level) %>%
  rowwise() %>%
  mutate(phe_center = round(mean(c(phe_from, phe_to)))) %>%
  ungroup()

cwindow_fdr_allchr <- list()
cpeak_fdr_allchr <- list()
for(chr in unique(cwindow_fdr$chr)){
  cwindow_fdr_chr <- filter(cwindow_fdr, phe_chr == !!chr)
  cpeak_fdr_chr <- filter(cpeak_fdr, phe_chr == !!chr)
  
  ## 2. add cols of if_overlap with the other phe & if_rep with significant signals in the other phe -----
  ## cwindow
  cwindow_fdr_allchr[[chr]] <- pbmcapply::pbmcmapply(
    cwindow_fdr_chr$phe_from,
    cwindow_fdr_chr$phe_to,
    cwindow_fdr_chr$phe_center,
    FUN = function(x, y, z){
      tmp_ind_overlap = data.table::between(pull(cpeak_fdr_chr, phe_center), x, y) |
        (
          (cpeak_fdr_chr$phe_from <= z) & (cpeak_fdr_chr$phe_to >= z)
        )
      
      tmp_ind_rep = cpeak_fdr_chr$if_sig & tmp_ind_overlap
      
      list(tibble(
        if_overlap = sum(tmp_ind_overlap),
        overlap_from = paste(cpeak_fdr_chr[tmp_ind_overlap, ]$phe_from, collapse = ';'),
        overlap_to = paste(cpeak_fdr_chr[tmp_ind_overlap, ]$phe_to, collapse = ';'),
        overlap_id = paste(cpeak_fdr_chr[tmp_ind_overlap, ]$phe_id, collapse = ';'),
        
        if_rep = sum(tmp_ind_rep),
        rep_from = paste(cpeak_fdr_chr[tmp_ind_rep, ]$phe_from, collapse = ';'),
        rep_to = paste(cpeak_fdr_chr[tmp_ind_rep, ]$phe_to, collapse = ';'),
        rep_id = paste(cpeak_fdr_chr[tmp_ind_rep, ]$phe_id, collapse = ';')
      ))
    },
    mc.cores = 10
  ) %>%
    bind_rows() %>%
    bind_cols(cwindow_fdr_chr)
  
  ## cpeak
  cpeak_fdr_allchr[[chr]] <- pbmcapply::pbmcmapply(
    cpeak_fdr_chr$phe_from,
    cpeak_fdr_chr$phe_to,
    cpeak_fdr_chr$phe_center,
    FUN = function(x, y, z){
      tmp_ind_overlap = data.table::between(pull(cwindow_fdr_chr, phe_center), x, y) |
        (
          (cwindow_fdr_chr$phe_from <= z) & (cwindow_fdr_chr$phe_to >= z)
        )
      
      tmp_ind_rep = cwindow_fdr_chr$if_sig & tmp_ind_overlap
      
      list(tibble(
        if_overlap = sum(tmp_ind_overlap),
        overlap_from = paste(cwindow_fdr_chr[tmp_ind_overlap, ]$phe_from, collapse = ';'),
        overlap_to = paste(cwindow_fdr_chr[tmp_ind_overlap, ]$phe_to, collapse = ';'),
        overlap_id = paste(cwindow_fdr_chr[tmp_ind_overlap, ]$group, collapse = ';'),
        
        if_rep = sum(tmp_ind_rep),
        rep_from = paste(cwindow_fdr_chr[tmp_ind_rep, ]$phe_from, collapse = ';'),
        rep_to = paste(cwindow_fdr_chr[tmp_ind_rep, ]$phe_to, collapse = ';'),
        rep_id = paste(cwindow_fdr_chr[tmp_ind_rep, ]$group, collapse = ';')
      ))
    },
    mc.cores = 10
  ) %>%
    bind_rows() %>%
    bind_cols(cpeak_fdr_chr)
}

cwindow_fdr_allchr <- bind_rows(cwindow_fdr_allchr)
cpeak_fdr_allchr <- bind_rows(cpeak_fdr_allchr)

# plot overlaps between cwindow and cpeak -----
# compute shared count once
n_shared <- cwindow_fdr_allchr %>%
  filter(if_rep > 0 & if_sig) %>%
  nrow()

plt_dat <- cwindow_fdr_allchr %>%
  mutate(
    method = ifelse(if_sig, "CACTI", "Single-peak based"),
    category = case_when(
      if_rep > 0 & if_sig  ~ "Shared",
      if_rep == 0 & if_sig ~ "Unique",
      if_rep > 0 & !if_sig ~ "Missed",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  dplyr::count(method, category) %>%
  add_row(
    method = "Single-peak based",
    category = "Shared",
    n = n_shared
  )

plt_dat$category <- factor(
  plt_dat$category,
  levels = c("Unique", "Missed", "Shared")
)

plt_dat$method <- factor(
  plt_dat$method,
  levels = c("CACTI", "Single-peak based")
)

ggplot(plt_dat, aes(x = method, y = n, fill = category)) +
  geom_col(width = 0.75, color = "black") +
  scale_fill_manual(
    values = c(
      Shared = "#5f9ea0",
      Unique = "#c6b7d9",
      Missed = "#f4b183"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of windows",
    fill = NULL
  ) +
  theme_classic(base_size = 7) +
  theme(
    legend.position = "right"
  )

