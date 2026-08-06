#!/usr/bin/r


# file: 2_classify_SNPs
# aim: classsify SNPs
# last update: 06-08-2026


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
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/")

snp_db <- fread("original_db/final_dataset_cna.csv")[, -1]
snp_metadata <- read_excel("original_db/Tabella pz TP53 (5).xlsx", sheet = "SNP")


# ---- data manipulation of the overall CNA dataset ----
# transposing the original dataframe
chr_arms <- mixedsort(unique(snp_db$chrarm))
  
snp_tmp_df <- snp_db %>% 
  select(N_SNP, SNP, MPC = Experiment_MPC_SNP.MPC, DNA, GITC, chrarm, weighted_mean_CN)

snp_tmp_df <- snp_tmp_df[mixedorder(snp_tmp_df$chrarm), ]
snp_tmp_df$chrarm <- factor(snp_tmp_df$chrarm, levels = unique(snp_tmp_df$chrarm))

snp_tmp_df <- snp_tmp_df %>% arrange(N_SNP)
snp_df <- dcast(setDT(snp_tmp_df), ... ~ chrarm, value.var = "weighted_mean_CN") %>% arrange(N_SNP)

# classify CN by arms - 10%
cna_10_df <- classifySNPArms(snp_df, 0.10, chr_arms)
write_tsv(cna_10_df, "data/cna/cna_class_10_df.txt")
write.xlsx(cna_10_df, "data/cna/cna_class_10_df.xlsx")

# classify CN by arms - 20%
cna_20_df <- classifySNPArms(snp_df, 0.20, chr_arms)
write_tsv(cna_20_df, "data/cna/cna_class_20_df.txt")
write.xlsx(cna_20_df, "data/cna/cna_class_20_df.xlsx")
    
# classify CN by arms - 50%
cna_50_df <- classifySNPArms(snp_df, 0.50, chr_arms)
write_tsv(cna_50_df, "data/cna/cna_class_50_df.txt")
write.xlsx(cna_50_df, "data/cna/cna_class_50_df.xlsx")

# classify CN by arms - 80%
cna_80_df <- classifySNPArms(snp_df, 0.80, chr_arms)
write_tsv(cna_80_df, "data/cna/cna_class_80_df.txt")
write.xlsx(cna_80_df, "data/cna/cna_class_80_df.xlsx")


# ---- data manipulation of filtered CNA dataset for 32 EMN02 pts ----
# homogenization of SNP ID
snp_meta <- snp_metadata %>% mutate(N_SNP = str_replace(SNP_numero, " ", "_"))
SNPs <- na.omit(snp_meta$N_SNP)

# 10% 
wb_cna_10_df <- cna_10_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_10_df, paste0("cna/wb_cna_10_df.txt"))
write.xlsx(wb_cna_10_df, "cna/wb_cna_10_df.xlsx")

# 20% 
wb_cna_20_df <- cna_20_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_20_df, paste0("cna/wb_cna_20_df.txt"))
write.xlsx(wb_cna_20_df, "cna/wb_cna_20_df.xlsx")

# 50% 
wb_cna_50_df <- cna_50_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_50_df, paste0("cna/wb_cna_50_df.txt"))
write.xlsx(wb_cna_50_df, "cna/wb_cna_50_df.xlsx")

# 80% 
wb_cna_80_df <- cna_80_df %>% 
  filter(N_SNP %in% SNPs)
write_tsv(wb_cna_80_df, paste0("cna/wb_cna_80_df.txt"))
write.xlsx(wb_cna_80_df, "cna/wb_cna_80_df.xlsx")


