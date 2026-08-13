# !/usr/bin/r


# file: 4_mmrf_isoform_combinations
# aim: check for co-expressed combinations of isoforms
# last update: 07-08-2026


library(openxlsx)
library(tidyverse)


# ---- Main ----
# setting environment
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/sanity_check/")


# ---- Cohort ----
mmrf_first_line_per_pt <- read.xlsx("mmrf_first_line_per_pt.xlsx")


# ---- Responses categories ----
mmrf_responses <- mmrf_first_line_per_pt %>%
  mutate(response_MU_PR = ifelse(resp_sh %in% c("PD", "SD"), 0, 1), 
         response_MU_VGPR = ifelse(resp_sh %in% c("PD", "SD", "PR"), 0, 1),
         response_MU_CR = ifelse(resp_sh %in% c("CR", "sCR"), 1, 0), 
         .after = resp_sh)


# ---- Isoforms combinations ----
mmrf_combinations <- mmrf_responses %>%
  mutate(p53_FL_1_0 = ifelse(p53_FL_exp=="absent", 0, 1), 
         p53_beta_1_0 = ifelse(p53_beta_exp=="absent", 0, 1),
         p53_gamma_1_0 = ifelse(p53_gamma_exp=="absent", 0, 1),
         delta40_p53_alpha_1_0 = ifelse(delta40_p53_alpha_exp=="absent", 0, 1),
         delta133_p53_alpha_1_0 = ifelse(delta133_p53_alpha_exp=="absent", 0, 1),
         delta133_p53_beta_1_0 = ifelse(delta133_p53_beta_exp=="absent", 0, 1),
         delta133_p53_gamma_1_0 = ifelse(delta133_p53_gamma_exp=="absent", 0, 1),
         FL_beta_133_alpha = ifelse(p53_FL_1_0==1 & p53_beta_1_0==1 & delta133_p53_alpha_1_0==1, 1, 0),
         FL_beta_gamma_133_alpha = ifelse(p53_FL_1_0==1 & p53_beta_1_0==1 & p53_gamma_1_0==1 & delta133_p53_alpha_1_0==1, 1, 0),
         FL_beta = ifelse(p53_FL_1_0==1 & p53_beta_1_0==1, 1, 0),
         FL_beta_gamma = ifelse(p53_FL_1_0==1 & p53_beta_1_0==1 & p53_gamma_1_0==1, 1, 0),
         beta_133_alpha = ifelse(p53_beta_1_0==1 & delta133_p53_alpha_1_0==1, 1, 0),
         .after = VAF) %>%
  arrange(PUBLIC_ID)

write.xlsx(mmrf_combinations, "9_isoform_combination_revision/mmrf_first_line_complete.xlsx")
