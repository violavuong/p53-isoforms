# !/usr/bin/r

# file: 2_mmrf_cna_categorization
# aim: categorize data using a 10% cut-off at two level, amp/del (1/0), and amp/del/normal
# last update: 13-08-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
# setting environment
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/genomics/")


# ---- Load Cohort ----
mmrf_gatk_cna_per_pt <- read.xlsx("mmrf_gatk_genome_per_pt_ia24.xlsx") %>%
  filter(line==1, specimen=="BM")


# ---- Categorization ----
mmrf_gatk_cna_categorized <- mmrf_gatk_cna_per_pt %>%
  mutate(across(c(4:48), ~ifelse(. > 2.10, 1, 0), .names = "{col}_amp"), 
         across(c(4:48), ~ifelse(. < 1.90, 1, 0), .names = "{col}_del"), 
         across(c(4:48), ~ifelse(. > 2.10, "amp", ifelse(. < 1.90, "del", "normal")), .names = "{col}_alt")) %>%
  select(-line, -specimen) %>%
  arrange(public_id)

mmrf_gatk_cna_categorized[is.na(mmrf_gatk_cna_categorized)] <- "NA"

write.xlsx(mmrf_gatk_cna_categorized, "mmrf_cna_per_pt.xlsx")
