##############################################
########### Enrichment of cqtls in eqtls ###########
##############################################
# load packages -----
rm(list = ls())
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

file_eqtlgen_best_hit <- '2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz'
fdr_level <- 0.05

# combine datasets
## LCL
mark_seq_lcl  <- 'H3K27AC|H3K4ME1|H3K4ME3'
list_file_lcl <- c(
  list.files(
    'results',
    str_glue('^cpeak_cwindow_({mark_seq_lcl})_window50kb_cis500kb.txt.gz$'),
    full.names = TRUE
  ),
  
  list.files(
    'results',
    str_glue('^cpeak_cwindow_seg_group_H3K36ME3_window50kb_cis500kb_seg5kb.txt.gz$'),
    full.names = TRUE
  )
)

## macrophage
mark_seq_macrophage  <- 'H3K27ac|H3K4me1|H3K4me3'
condition_macrophage  <- 'NI|Flu'
list_file_macrophage <- c(
  list.files(
    'results',
    str_glue('^cpeak_cwindow_({mark_seq_macrophage})_({condition_macrophage})_window50kb_cis100kb.txt.gz$'),
    full.names = TRUE
  ),
  
  list.files(
    'results',
    str_glue('^cpeak_cwindow_seg_group_H3K27me3_({condition_macrophage})_window50kb_seg5kb_cis100kb.txt.gz$'),
    full.names = TRUE
  )
)

## output -----
lcl_outputs <- c(
  str_glue('{list_file_lcl}_new_cqtl_eqtl_enrich.pdf'),
  str_glue('{list_file_lcl}_all_cqtl_eqtl_enrich.pdf'),
  str_glue('{list_file_lcl}_shared_cqtl_eqtl_enrich.pdf')
)
macrophage_outputs <- c(
  str_glue('{list_file_macrophage}_new_cqtl_eqtl_enrich.pdf'),
  str_glue('{list_file_macrophage}_all_cqtl_eqtl_enrich.pdf'),
  str_glue('{list_file_macrophage}_shared_cqtl_eqtl_enrich.pdf')
)
combined_output <- c(
  'cqtl_eqtl_enrich_OR_across_datasets.pdf',
  'cqtl_new_eqtl_enrich_across_datasets.pdf',
  'cqtl_all_eqtl_enrich_across_datasets.pdf'
)

all_outputs <- c(lcl_outputs, macrophage_outputs, combined_output)

## Print ----
cat(
  "=== Figure outputs that will be generated ===\n",
  "Total figures:", length(all_outputs), "\n"
)
cat(all_outputs, sep = " ")



# read files -----
## eqtlgen -----
eqtlgen_best_hit <- data.table::fread(file_eqtlgen_best_hit)


# extra new cqtls overlap with eqtlgen lead esnps -----
## take eqtlgen eqtls -----
## all eqtlgen eqtls
eqtlgen_best_hit$snp_hg19 <- paste(eqtlgen_best_hit$SNPChr, eqtlgen_best_hit$SNPPos, sep = ":")

## only lead snp for each gene
eqtlgen_top <- group_by(eqtlgen_best_hit, GeneSymbol) %>%
  summarise(
    best_hit = min(Pvalue),
    n_var_in_cis = n(),
    snp_hg19 = snp_hg19[which.min(Pvalue)]
  ) %>%
  ungroup()


## lcl -----
dict <- c(
  H3K27AC  = "H3K27ac",
  H3K4ME1  = "H3K4me1",
  H3K4ME3  = "H3K4me3",
  H3K36ME3 = "H3K36me3"
)

window_best_hit_lcl <- lapply(
  list_file_lcl,
  function(x) {
    mark = dict[str_extract(x, 'H3K27AC|H3K4ME1|H3K4ME3|H3K36ME3')]
    
    window_best_hit_tmp = data.table::fread(x)
    
    ## output -----
    file_bed_window_cqtl = str_glue('{x}_cqtl.bed')
    file_conversion = str_glue('{x}_conversion.bed')
    file_unmap = str_glue('{x}_unmap')
    
    
    # liftover to hg19 -----
    ## match SNPs -----
    conversion = data.table::fread(
      file_conversion, 
      header = FALSE, 
      col.names = c('chr', 'pos_l', 'pos_hg19', 'snp')
    )
    
    
    ## change line of snp_hg19, mutate conversion to snp_hg19 first, in case conversion file changes chr
    window_best_hit_tmp = left_join(
      window_best_hit_tmp, select(conversion, snp, pos_hg19),
      by = c('snp')
    ) %>%
      mutate(
        snp_hg19 = case_when(!is.na(pos_hg19) ~ paste(str_extract(chr, '\\d+'), pos_hg19, sep = ":"), .default = NA),
        if_match = !is.na(pos_hg19)
      )
    
    
    ## plot single dataset
    ## simple overlap
    window_best_hit_tmp$if_eqtl_overlap <- window_best_hit_tmp$snp_hg19 %in% eqtlgen_best_hit$snp_hg19
    window_best_hit_tmp$if_eqtl_top <- window_best_hit_tmp$snp_hg19 %in% eqtlgen_top$snp_hg19
    
    # enrichment of extra new cqtls (combine all datasets) -----
    # new cqtls
    dat_extra_enrich_tmp_new_cqtl = bind_rows(
      filter(
        window_best_hit_tmp,
        if_match & if_sig & (if_rep == 0)
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_extra = TRUE
        ),
      
      filter(
        window_best_hit_tmp,
        if_match & !if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_extra = FALSE
        )
    ) %>%
      arrange(desc(if_extra), desc(if_eqtl_overlap))
    dat_extra_enrich_tmp_new_cqtl$if_extra <- ifelse(dat_extra_enrich_tmp_new_cqtl$if_extra, "New cQTL", "Not cQTL")
    dat_extra_enrich_tmp_new_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_tmp_new_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")
    
    # all cqtls
    dat_extra_enrich_tmp_all_cqtl = bind_rows(
      filter(
        window_best_hit_tmp,
        if_match & if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_cqtl = TRUE
        ),
      
      filter(
        window_best_hit_tmp,
        if_match & !if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_cqtl = FALSE
        )
    ) %>%
      arrange(desc(if_cqtl), desc(if_eqtl_overlap))
    dat_extra_enrich_tmp_all_cqtl$if_cqtl <- ifelse(dat_extra_enrich_tmp_all_cqtl$if_cqtl, "cQTL", "Not cQTL")
    dat_extra_enrich_tmp_all_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_tmp_all_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")
    
    
    # shared cqtls
    dat_extra_enrich_tmp_shared_cqtl = bind_rows(
      filter(
        window_best_hit_tmp,
        if_match & if_sig & (if_rep > 0)
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_shared = TRUE
        ),
      
      filter(
        window_best_hit_tmp,
        if_match & !if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_shared = FALSE
        )
    ) %>%
      arrange(desc(if_shared), desc(if_eqtl_overlap))
    dat_extra_enrich_tmp_shared_cqtl$if_shared <- ifelse(dat_extra_enrich_tmp_shared_cqtl$if_shared, "Shared cQTL", "Not cQTL")
    dat_extra_enrich_tmp_shared_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_tmp_shared_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")
    
    
    # test -----
    enrich_p_tmp_new_cqtl = fisher.test(matrix(dat_extra_enrich_tmp_new_cqtl$n, ncol = 2), alternative = "greater")
    enrich_p_tmp_all_cqtl = fisher.test(matrix(dat_extra_enrich_tmp_all_cqtl$n, ncol = 2), alternative = "greater")
    enrich_p_tmp_shared_cqtl = fisher.test(matrix(dat_extra_enrich_tmp_shared_cqtl$n, ncol = 2), alternative = "greater")
    
    
    ## plot enrichment -----
    # new cqtls
    ggplot(dat_extra_enrich_tmp_new_cqtl, aes(x = if_extra, y = n, fill = if_eqtl_overlap)) +
      geom_col() +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
      labs(
        x = NULL, y = "Number of cWindows", fill = NULL,
        title = str_glue(
          "Enrichment of new cQTLs for eQTLs\n {mark} in LCL \n P:{sprintf('%.2e', enrich_p_tmp_new_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_tmp_new_cqtl$estimate)}"
        )
      ) +
      scale_fill_manual(
        values = c(
          "eQTL"     = "#6bbf8a",
          "Not eQTL" = "#DDDDDD"
        )
      )
    
    ## save figure -----
    ggsave(
      str_glue('{x}_new_cqtl_eqtl_enrich.pdf'),
      width = 3, height = 2
    )
    
    
    ## plot enrichment -----
    # all cqtls
    ggplot(dat_extra_enrich_tmp_all_cqtl, aes(x = if_cqtl, y = n, fill = if_eqtl_overlap)) +
      geom_col() +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
      labs(
        x = NULL, y = "Number of cWindows", fill = NULL,
        title = str_glue(
          "Enrichment of all cQTLs for eQTLs\n {mark} in LCL \n P:{sprintf('%.2e', enrich_p_tmp_all_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_tmp_all_cqtl$estimate)}"
        )
      ) +
      scale_fill_manual(
        values = c(
          "eQTL"     = "#6bbf8a",
          "Not eQTL" = "#DDDDDD"
        )
      )
    
    ## save figure -----
    ggsave(
      str_glue('{x}_all_cqtl_eqtl_enrich.pdf'),
      width = 3, height = 2
    )
    
    # shared cqtls
    ggplot(dat_extra_enrich_tmp_shared_cqtl, aes(x = if_shared, y = n, fill = if_eqtl_overlap)) +
      geom_col() +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
      labs(
        x = NULL, y = "Number of cWindows", fill = NULL,
        title = str_glue(
          "Enrichment of shared cQTLs for eQTLs\n {mark} in LCL \n P:{sprintf('%.2e', enrich_p_tmp_shared_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_tmp_shared_cqtl$estimate)}"
        )
      ) +
      scale_fill_manual(
        values = c(
          "eQTL"     = "#6bbf8a",
          "Not eQTL" = "#DDDDDD"
        )
      )
    
    ## save figure -----
    ggsave(
      str_glue('{x}_shared_cqtl_eqtl_enrich.pdf'),
      width = 3, height = 2
    )
    
    
    
    return(list(
      'window_best_hit_tmp' = window_best_hit_tmp,
      'mark' = mark,
      'cell type' = 'LCL',
      'P (new cQTLs)' = enrich_p_tmp_new_cqtl$p.value, 'OR (new cQTLs)' = enrich_p_tmp_new_cqtl$estimate,
      'P (all cQTLs)' = enrich_p_tmp_all_cqtl$p.value, 'OR (all cQTLs)' = enrich_p_tmp_all_cqtl$estimate,
      'P (shared cQTLs)' = enrich_p_tmp_shared_cqtl$p.value, 'OR (shared cQTLs)' = enrich_p_tmp_shared_cqtl$estimate
    ))
  }
)


## macrophage -----
window_best_hit_macrophage <- lapply(
  list_file_macrophage,
  function(x) {
    mark = str_extract(x, '(H3K27ac|H3K4me1|H3K4me3|H3K27me3)_(Flu|NI)')
    
    window_best_hit_tmp = data.table::fread(x)
    window_best_hit_tmp$if_match = TRUE
    window_best_hit_tmp$snp_hg19 = str_extract_all(window_best_hit_tmp$snp, "\\d+") %>% sapply(paste, collapse = ":")
    
    ## plot single dataset
    ## simple overlap
    window_best_hit_tmp$if_eqtl_overlap <- window_best_hit_tmp$snp_hg19 %in% eqtlgen_best_hit$snp_hg19
    window_best_hit_tmp$if_eqtl_top <- window_best_hit_tmp$snp_hg19 %in% eqtlgen_top$snp_hg19
    
    # enrichment of extra new cqtls (combine all datasets) -----
    ## here use all eqtlgen eqtls and simple overlap
    # new cqtls
    dat_extra_enrich_tmp_new_cqtl = bind_rows(
      filter(
        window_best_hit_tmp,
        if_match & if_sig & (if_rep == 0)
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_extra = TRUE
        ),
      
      filter(
        window_best_hit_tmp,
        if_match & !if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_extra = FALSE
        )
    ) %>%
      arrange(desc(if_extra), desc(if_eqtl_overlap))
    dat_extra_enrich_tmp_new_cqtl$if_extra <- ifelse(dat_extra_enrich_tmp_new_cqtl$if_extra, "New cQTL", "Not cQTL")
    dat_extra_enrich_tmp_new_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_tmp_new_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")
    
    # all cqtls
    dat_extra_enrich_tmp_all_cqtl = bind_rows(
      filter(
        window_best_hit_tmp,
        if_match & if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_cqtl = TRUE
        ),
      
      filter(
        window_best_hit_tmp,
        if_match & !if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_cqtl = FALSE
        )
    ) %>%
      arrange(desc(if_cqtl), desc(if_eqtl_overlap))
    dat_extra_enrich_tmp_all_cqtl$if_cqtl <- ifelse(dat_extra_enrich_tmp_all_cqtl$if_cqtl, "cQTL", "Not cQTL")
    dat_extra_enrich_tmp_all_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_tmp_all_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")
    
    
    # shared cqtls
    dat_extra_enrich_tmp_shared_cqtl = bind_rows(
      filter(
        window_best_hit_tmp,
        if_match & if_sig & (if_rep > 0)
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_shared = TRUE
        ),
      
      filter(
        window_best_hit_tmp,
        if_match & !if_sig
      ) %>%
        count(if_eqtl_overlap) %>%
        mutate(
          if_shared = FALSE
        )
    ) %>%
      arrange(desc(if_shared), desc(if_eqtl_overlap))
    dat_extra_enrich_tmp_shared_cqtl$if_shared <- ifelse(dat_extra_enrich_tmp_shared_cqtl$if_shared, "Shared cQTL", "Not cQTL")
    dat_extra_enrich_tmp_shared_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_tmp_shared_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")
    
    
    # test -----
    enrich_p_tmp_new_cqtl = fisher.test(matrix(dat_extra_enrich_tmp_new_cqtl$n, ncol = 2), alternative = "greater")
    enrich_p_tmp_all_cqtl = fisher.test(matrix(dat_extra_enrich_tmp_all_cqtl$n, ncol = 2), alternative = "greater")
    enrich_p_tmp_shared_cqtl = fisher.test(matrix(dat_extra_enrich_tmp_shared_cqtl$n, ncol = 2), alternative = "greater")
    
    
    
    ## plot enrichment -----
    # new cqtls
    ggplot(dat_extra_enrich_tmp_new_cqtl, aes(x = if_extra, y = n, fill = if_eqtl_overlap)) +
      geom_col() +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
      labs(
        x = NULL, y = "Number of cWindows", fill = NULL,
        title = str_glue(
          "Enrichment of new cQTLs for eQTLs\n {mark} in Macrophage \n P:{sprintf('%.2e', enrich_p_tmp_new_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_tmp_new_cqtl$estimate)}"
        )
      ) +
      scale_fill_manual(
        values = c(
          "eQTL"     = "#6bbf8a",
          "Not eQTL" = "#DDDDDD"
        )
      )
    
    ## save figure -----
    ggsave(
      str_glue('{x}_new_cqtl_eqtl_enrich.pdf'),
      width = 3, height = 2
    )
    
    
    ## plot enrichment -----
    # all cqtls
    ggplot(dat_extra_enrich_tmp_all_cqtl, aes(x = if_cqtl, y = n, fill = if_eqtl_overlap)) +
      geom_col() +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
      labs(
        x = NULL, y = "Number of cWindows", fill = NULL,
        title = str_glue(
          "Enrichment of all cQTLs for eQTLs\n {mark} in Macrophage \n P:{sprintf('%.2e', enrich_p_tmp_all_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_tmp_all_cqtl$estimate)}"
        )
      ) +
      scale_fill_manual(
        values = c(
          "eQTL"     = "#6bbf8a",
          "Not eQTL" = "#DDDDDD"
        )
      )
    
    ## save figure -----
    ggsave(
      str_glue('{x}_all_cqtl_eqtl_enrich.pdf'),
      width = 3, height = 2
    )
    
    
    # shared cqtls
    ggplot(dat_extra_enrich_tmp_shared_cqtl, aes(x = if_shared, y = n, fill = if_eqtl_overlap)) +
      geom_col() +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
      labs(
        x = NULL, y = "Number of cWindows", fill = NULL,
        title = str_glue(
          "Enrichment of shared cQTLs for eQTLs\n {mark} in Macrophage \n P:{sprintf('%.2e', enrich_p_tmp_shared_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_tmp_shared_cqtl$estimate)}"
        )
      ) +
      scale_fill_manual(
        values = c(
          "eQTL"     = "#6bbf8a",
          "Not eQTL" = "#DDDDDD"
        )
      )
    
    ## save figure -----
    ggsave(
      str_glue('{x}_shared_cqtl_eqtl_enrich.pdf'),
      width = 3, height = 2
    )
    
    
    return(list(
      'window_best_hit_tmp' = window_best_hit_tmp,
      'mark' = mark,
      'cell type' = 'macrophage',
      'P (new cQTLs)' = enrich_p_tmp_new_cqtl$p.value, 'OR (new cQTLs)' = enrich_p_tmp_new_cqtl$estimate,
      'P (all cQTLs)' = enrich_p_tmp_all_cqtl$p.value, 'OR (all cQTLs)' = enrich_p_tmp_all_cqtl$estimate,
      'P (shared cQTLs)' = enrich_p_tmp_shared_cqtl$p.value, 'OR (shared cQTLs)' = enrich_p_tmp_shared_cqtl$estimate
    ))
  }
)


## enrich p and OR for single datasets -----
enrich_num_comb <- bind_rows(
  lapply(window_best_hit_lcl, FUN = function(x) x[names(x) != 'window_best_hit_tmp']),
  lapply(window_best_hit_macrophage, FUN = function(x) x[names(x) != 'window_best_hit_tmp'])
) %>%
  mutate(
    across(starts_with("P ("),  ~ sprintf("%.2e", as.numeric(.x))),
    across(starts_with("OR ("), ~ sprintf("%.2f", as.numeric(.x)))
  )

### vis -----
# Define significance thresholds
# *: P < 0.05, **: P < 0.01, ***: P < 0.001
enrich_num_comb %>%
  # 1. Pivot OR and P columns to long format first
  pivot_longer(
    cols = matches("^(OR|P)"), 
    names_to = "raw_name", 
    values_to = "value"
  ) %>%
  # 2. Extract the Metric (OR vs P) and the Type (all/shared/new)
  mutate(
    Metric = str_extract(raw_name, "^(OR|P)"),
    Type = str_extract(raw_name, "(?<= \\().*(?=\\))"), # Extracts text inside parenthesis
    value = as.numeric(value) # Ensure numbers are numeric
  ) %>%
  select(-raw_name) %>%
  # 3. Spread OR and P back into separate columns so they are side-by-side
  pivot_wider(
    names_from = Metric,
    values_from = value
  ) %>%
  # 4. Add Stars
  mutate(
    Type = factor(Type, levels = c("all cQTLs", "shared cQTLs", "new cQTLs")),
    stars = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ "" # No star if not significant
    )
  ) %>%
  # 5. Plot
  ggplot(aes(x = str_glue('{mark}\n({`cell type`})'), y = OR, fill = Type)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(
    aes(label = stars),
    position = position_dodge2(width = 0.9, preserve = "single"),
    vjust = -0.2,
    size = 8 / .pt, 
    fontface = "bold"
  ) +
labs(
  x = "Dataset", y = "Odds Ratio (OR)", fill = "cQTL type",
  title = "Enrichment of cQTLs for eQTLs across datasets"
) +
  scale_fill_manual(
    values = c(
      "all cQTLs"    = "#e49e61",
      "shared cQTLs" = "#4c9b9b",
      "new cQTLs"    = "#9575ab"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


## combine all datasets -----
window_best_hit <- bind_rows(
  lapply(window_best_hit_lcl, FUN = function(x) x$window_best_hit_tmp),
  lapply(window_best_hit_macrophage, FUN = function(x) x$window_best_hit_tmp)
) %>%
  group_by(snp_hg19) %>%
  summarise(
    if_match = any(if_match),
    if_sig = any(if_sig),
    if_rep = any(if_rep)
  ) %>%
  ungroup()


## Define overlap -----
## simple overlap
window_best_hit$if_eqtl_overlap <- window_best_hit$snp_hg19 %in% eqtlgen_best_hit$snp_hg19
window_best_hit$if_eqtl_top <- window_best_hit$snp_hg19 %in% eqtlgen_top$snp_hg19



# enrichment of cqtls (combine all datasets) -----
# new cqtls
dat_extra_enrich_new_cqtl <- bind_rows(
  filter(
    window_best_hit,
    if_match & if_sig & (if_rep == 0)
  ) %>%
    count(if_eqtl_overlap) %>%
    mutate(
      if_extra = TRUE
    ),
  
  filter(
    window_best_hit,
    if_match & !if_sig
  ) %>%
    count(if_eqtl_overlap) %>%
    mutate(
      if_extra = FALSE
    )
) %>%
  arrange(desc(if_extra), desc(if_eqtl_overlap))
dat_extra_enrich_new_cqtl$if_extra <- ifelse(dat_extra_enrich_new_cqtl$if_extra, "New cQTL", "Not cQTL")
dat_extra_enrich_new_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_new_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")

## test -----
enrich_p_new_cqtl <- fisher.test(matrix(dat_extra_enrich_new_cqtl$n, ncol = 2), alternative = "greater")


## plot enrichment -----
ggplot(dat_extra_enrich_new_cqtl, aes(x = if_extra, y = n, fill = if_eqtl_overlap)) +
  geom_col() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
  labs(
    x = NULL, y = "Number of cWindows", fill = NULL,
    title = str_glue(
      "Enrichment of new cQTLs for eQTLs\n All marks \n P:{sprintf('%.2e', enrich_p_new_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_new_cqtl$estimate)}"
    )
  ) +
  scale_fill_manual(
    values = c(
      "eQTL"     = "#6bbf8a",
      "Not eQTL" = "#DDDDDD"
    )
  )


# all cqtls
dat_extra_enrich_all_cqtl = bind_rows(
  filter(
    window_best_hit,
    if_match & if_sig
  ) %>%
    count(if_eqtl_overlap) %>%
    mutate(
      if_cqtl = TRUE
    ),
  
  filter(
    window_best_hit,
    if_match & !if_sig
  ) %>%
    count(if_eqtl_overlap) %>%
    mutate(
      if_cqtl = FALSE
    )
) %>%
  arrange(desc(if_cqtl), desc(if_eqtl_overlap))
dat_extra_enrich_all_cqtl$if_cqtl <- ifelse(dat_extra_enrich_all_cqtl$if_cqtl, "cQTL", "Not cQTL")
dat_extra_enrich_all_cqtl$if_eqtl_overlap <- ifelse(dat_extra_enrich_all_cqtl$if_eqtl_overlap, "eQTL", "Not eQTL")

## test -----
enrich_p_all_cqtl <- fisher.test(matrix(dat_extra_enrich_all_cqtl$n, ncol = 2), alternative = "greater")


## plot enrichment -----
# all cqtls
ggplot(dat_extra_enrich_all_cqtl, aes(x = if_cqtl, y = n, fill = if_eqtl_overlap)) +
  geom_col() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6 / .pt) +
  labs(
    x = NULL, y = "Number of cWindows", fill = NULL,
    title = str_glue(
      "Enrichment of all cQTLs for eQTLs\n All marks \n P:{sprintf('%.2e', enrich_p_all_cqtl$p.value)}, OR:{sprintf('%.2f', enrich_p_all_cqtl$estimate)}"
    )
  ) +
  scale_fill_manual(
    values = c(
      "eQTL"     = "#6bbf8a",
      "Not eQTL" = "#DDDDDD"
    )
  )

