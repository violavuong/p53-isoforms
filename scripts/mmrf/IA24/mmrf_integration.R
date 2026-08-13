# !/usr/bin/r


# file: mmrf_integration
# aim: merge genomic and transcriptomic data based on the key "public_id"
# last update: 13-08-26


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
