# !/usr/bin/r


# file: 1_check_datasets
# aim: checking data consistency
# last update: 06-08-2026


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
# environment setting
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")
outDir <- "/mmrf/sanity_check/"

dir.create(paste0(wd, outDir), recursive = TRUE)


# ---- Sanity check: transcript-based ----
# loading original mmrf dataset
mmrf_transcript_based_df <- fread(paste0(wd, "/transcript_based/mmrf_tp53_per_pt.txt"))
mmrf_genomic_based_df <- fread(paste0(wd, "/genomic_based/mmrf_tp53_per_pt.txt"))