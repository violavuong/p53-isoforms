# !/usr/bin/r


# file: 5_mmrf_second_line_survival
# aim: check ttp, and PFS events of pts who progressed to second line
# last update: 10-08-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
# setting environment
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/sanity_check/")
outDir <- "/survival_second_line/"

dir.create(paste0(wd, outDir), recursive = TRUE)


# ---- Cohort ----
mmrf_complete_df <- read.xlsx("mmrf_first_line_per_pt_complete.xlsx")


# ---- Check survival data of MMRF IA22 ----
mmrf_survival <- fread("survival_second_line/MMRF_CoMMpass_IA22_STAND_ALONE_SURVIVAL.tsv") %>%
  filter(PUBLIC_ID %in% mmrf_complete_df$PUBLIC_ID) %>%
  arrange(PUBLIC_ID)

write.xlsx(mmrf_survival, "mmrf_survival_659_pts.xlsx")

# ---- Merge data ----
mmrf_complete_with_ttp_df <- left_join(mmrf_complete_df, mmrf_survival_second_line, by = "PUBLIC_ID") %>%
  relocate(c(start_line_2, end_line_2, best_resp_dy_line_2), .after = best_resp_dy_line_1) %>%
  relocate(c(PFS_date_line_2, PFS_event_line_2, PFS_time_line_2, ttp, ttp_event), .after = PFS_months) %>%
  arrange(PUBLIC_ID)

write.xlsx(mmrf_complete_with_ttp_df, "mmrf_first_line_per_pt_with_ttp.xlsx")  
