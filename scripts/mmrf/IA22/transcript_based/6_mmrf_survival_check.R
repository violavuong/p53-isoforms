# !/usr/bin/r


# file: 6_mmrf_survival_check
# aim: attach correct survival columns
# last update: 10-08-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
# setting environment 
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/sanity_check/")


# ---- Cohort ----
mmrf_complete_df <- read.xlsx("mmrf_first_line_per_pt_complete.xlsx") %>%
  select(-best_resp_dy_line_1, -starts_with("PFS"), -starts_with("OS"))


# ---- Survival data ----
mmrf_survival_per_pt <- fread("survival_second_line/MMRF_CoMMpass_IA22_STAND_ALONE_SURVIVAL.tsv") %>%
  filter(PUBLIC_ID %in% mmrf_complete_df$PUBLIC_ID) %>%
  select(PUBLIC_ID, deathdy, lstalive, pdflag, ttfpd, censpfs, pfscdy, censos, oscdy, ttcos)


# ---- Merge ----
mmrf_correct_survival_per_pt <- left_join(mmrf_complete_df, mmrf_survival_per_pt, by = "PUBLIC_ID") %>%
  relocate(c(deathdy, lstalive, pdflag, ttfpd, censpfs, pfscdy, censos, oscdy, ttcos), .after = end_line_1) %>%
  arrange(PUBLIC_ID)

write.xlsx(mmrf_correct_survival_per_pt, "mmrf_first_line_per_pt_corrected.xlsx")
