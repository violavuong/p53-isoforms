#!/usr/bin/r

## file: SNP_classification
## last update: 20-06-2025

# ---- Main ----
library(data.table)
library(gtools)
library(readxl)
library(purrr)
library(tidyverse)


# global input
setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/")

snp_db <- fread("data/original_db/final_dataset_cna.csv")[, -1]
snp_meta <- read_excel("data/original_db/Tabella pz TP53 (5).xlsx", sheet = "SNP")

# ---- data manipulation of overall CNA dataset ----
snp_db <- snp_db %>% mutate(abs_CN_diff = abs(weighted_mean_CN - 2))

cut_offs_cols <- map2_dfc(.x = names(cut_offs), 
                          .y = as.list(cut_offs), 
                          .f = ~ {tibble(!!sym(.x) := if_else(snp_db$abs_CN_diff > .y, 1, 0))})

snp_db <- cbind(snp_db, cut_offs_cols)
snp_db <- snp_db %>% relocate(abs_CN_diff, class_10, class_20, class_50, class_80, .after = weighted_mean_CN)
write_tsv(snp_db, paste0("data/cna_class_df.txt"))


# ---- data manipulation of filtered CNA dataset for 32 EMN02 pts ----
#homogenization of SNP ID
snp_meta <- snp_meta %>% mutate(N_SNP = str_replace(SNP_numero, " ", "_"))
SNPs <- na.omit(snp_meta$N_SNP)

#filtering and arranging by chrarm and SNP ID
wb_snp_tmp_df <- snp_db %>% 
  filter(N_SNP %in% SNPs) %>%
  select(MPC = Experiment_MPC_SNP.MPC, DNA, GITC, N_SNP, SNP, chrarm, weighted_mean_CN)
wb_snp_tmp_df <- wb_snp_tmp_df[mixedorder(wb_snp_tmp_df$chrarm), ]
wb_snp_tmp_df <- wb_snp_tmp_df %>% arrange(N_SNP)


# ---- classifying by selected cut-off ----
cut_offs <- c(class_10 = 0.10, class_20 = 0.20, class_50 = 0.50, class_80 = 0.80)
wb_snp_tmp_df <- wb_snp_tmp_df %>% mutate(abs_CN_diff = abs(weighted_mean_CN - 2))

cut_offs_cols <- map2_dfc(.x = names(cut_offs), 
                          .y = as.list(cut_offs), 
                          .f = ~ {tibble(!!sym(.x) := if_else(wb_snp_tmp_df$abs_CN_diff > .y, 1, 0))})

wb_snp_df <- cbind(wb_snp_tmp_df, cut_offs_cols)

write_tsv(wb_snp_df, paste0("data/wb_snp_df.txt"))

