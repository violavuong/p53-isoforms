#!/usr/bin/r

## file: SNP_classification
## last update: 25-06-2025

# ---- Main ----
library(data.table)
library(gtools)
library(openxlsx)
library(readxl)
library(rlist)
library(purrr)
library(tidyverse)

source("script/fun/utils.R")

# global input
setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/")

cna_db <- fread("data/original_db/final_dataset_cna.csv")[, -1]
snp_meta <- read_excel("data/original_db/Tabella pz TP53 (5).xlsx", sheet = "SNP")


# ---- data manipulation of the overall CNA dataset ----
# transposing the original dataframe
chr_arms <- mixedsort(unique(cna_db$chrarm))
  
cna_tmp_df <- cna_db %>% 
  select(N_SNP, SNP, MPC = Experiment_MPC_SNP.MPC, DNA, GITC, chrarm, weighted_mean_CN)
cna_tmp_df <- cna_tmp_df[mixedorder(cna_tmp_df$chrarm), ]
cna_tmp_df$chrarm <- factor(cna_tmp_df$chrarm, levels = unique(cna_tmp_df$chrarm))
cna_tmp_df <- cna_tmp_df %>% arrange(N_SNP)
cna_df <- dcast(setDT(cna_tmp_df), ... ~ chrarm, value.var = "weighted_mean_CN") %>% arrange(N_SNP)

# classify CN by arms - 10%
cna_10_df <- classifyArms(cna_df, 0.1, chr_arms)
write_tsv(cna_10_df, "data/cna/cna_class_10_df.txt")
write.xlsx(cna_10_df, "data/cna/cna_class_10_df.xlsx")

# classify CN by arms - 20%
cna_20_df <- classifyArms(cna_df, 0.2, chr_arms)
write_tsv(cna_20_df, "data/cna/cna_class_20_df.txt")
write.xlsx(cna_20_df, "data/cna/cna_class_20_df.xlsx")
    
# classify CN by arms - 50%
cna_50_df <- classifyArms(cna_df, 0.5, chr_arms)
write_tsv(cna_50_df, "data/cna/cna_class_50_df.txt")
write.xlsx(cna_50_df, "data/cna/cna_class_50_df.xlsx")

# classify CN by arms - 80%
cna_80_df <- classifyArms(cna_df, 0.8, chr_arms)
write_tsv(cna_80_df, "data/cna/cna_class_80_df.txt")
write.xlsx(cna_80_df, "data/cna/cna_class_80_df.xlsx")


# ---- data manipulation of filtered CNA dataset for 32 EMN02 pts ----
# homogenization of SNP ID
snp_meta <- snp_meta %>% mutate(N_SNP = str_replace(SNP_numero, " ", "_"))
SNPs <- na.omit(snp_meta$N_SNP)

# 10% 
wb_cna_10_df <- cna_10_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_10_df, paste0("data/cna/wb_cna_10_df.txt"))
write.xlsx(wb_cna_10_df, "data/cna/wb_cna_10_df.xlsx")

# 20% 
wb_cna_20_df <- cna_20_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_20_df, paste0("data/cna/wb_cna_20_df.txt"))
write.xlsx(wb_cna_20_df, "data/cna/wb_cna_20_df.xlsx")

# 50% 
wb_cna_50_df <- cna_50_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_50_df, paste0("data/cna/wb_cna_50_df.txt"))
write.xlsx(wb_cna_50_df, "data/cna/wb_cna_50_df.xlsx")

# 80% 
wb_cna_80_df <- cna_80_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_80_df, paste0("data/cna/wb_cna_80_df.txt"))
write.xlsx(wb_cna_80_df, "data/cna/wb_cna_80_df.xlsx")






