# !/usr/bin/r


# file: mmrf_integration
# aim: merge genomic and transcriptomic data based on the key "public_id"
# last update: 02-09-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/")


# ---- Loading datasets ----
mmrf_tp53_per_pt <- read.xlsx("transcriptomic/mmrf_tp53_per_pt.xlsx")
  
mmrf_cna_per_pt <- read.xlsx("genomics/mmrf_cna_per_pt.xlsx")

mmrf_translocations_per_pt <- read.xlsx("genomics/mmrf_translocation_per_pt.xlsx")

mmrf_mut_per_pt <- fread("genomics/mmrf_ns_mut_per_pt.txt") %>%
  select(public_id, TP53_mut = TP53)


# ---- Merge ----
mmrf_per_pt <- inner_join(mmrf_tp53_per_pt, mmrf_cna_per_pt, by = "public_id") %>%
  left_join(mmrf_translocations_per_pt, by = "public_id") %>%
  left_join(mmrf_mut_per_pt, by = "public_id") %>%
  arrange(public_id)

mmrf_per_pt[is.na(mmrf_per_pt)] <- "na"

write.xlsx(mmrf_per_pt, "mmrf_omics_per_pt.xlsx")


# ---- adding clinical data ----
mmrf_per_pt <- read.xlsx("mmrf_omics_per_pt.xlsx")

mmrf_surv <- read.xlsx("commpass_surv.xlsx")

mmrf_biochem <- read.xlsx("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/MM group - Vuong_Viola_Meixian/ims_analysis/old_mmrf_IA24/db_analisi/clinical/mmrf_biochem_per_pt_030626.xlsx")

mmrf_labs <- fread("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/clinical_data/labs_deid.csv")

mmrf_complete_per_pt <- left_join(mmrf_per_pt, mmrf_surv, by = "public_id") %>%
  left_join(mmrf_biochem, by = "public_id")

write.xlsx(mmrf_complete_per_pt, "mmrf_complete_3_data_per_pt.xlsx")
