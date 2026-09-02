# !/usr/bin/r


# file: 2_mmrf_cohort_definition
# aim: filtering the original dataset for first line pts expressing selected p53 isoforms
# last update: 02-09-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/transcriptomic/")
mmrfDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/transcript_expression/"


# ---- Define cohort ----
mmrf_p53_tpm_per_pt <- fread("complete_db/mmrf_tpm_per_pt_ia24.txt") %>%
  filter(line==1, specimen=="BM") %>%
  select(public_id, p53_FL = ENST00000269305, p53_beta = ENST00000420246, p53_gamma = ENST00000455263, delta40_p53_alpha = ENST00000610292, 
         delta133_p53_alpha = ENST00000504937, delta133_p53_beta = ENST00000510385, delta133_p53_gamma = ENST00000504290)


# ---- Categorization ----
# binary categorization 
mmrf_p53_presence_cat <- mmrf_p53_tpm_per_pt %>%
  mutate(across(c(2:8), ~ifelse(.==0, 0, 1), .names = "{col}_1_0"))

# expression categorization 
median_per_isoform <- mmrf_p53_tpm_per_pt %>% 
  summarise(across(c(2:8), ~median(.[. > 0])))

mmrf_p53_expression_cat <- mmrf_p53_presence_cat %>%
  mutate(across(c(2:8), ~ifelse(.==0, "absent", ifelse(.x > median_per_isoform[[cur_column()]], "high", "low")), .names = "{col}_exp")) %>%
  select(public_id, starts_with("p53_FL"), starts_with("p53_beta"), starts_with("p53_gamma"), starts_with("delta40"), 
         starts_with("delta133_p53_alpha"), starts_with("delta133_p53_beta"), starts_with("delta133_p53_gamma")) %>%
  arrange(public_id)


write.xlsx(mmrf_p53_expression_cat, "mmrf_tp53_per_pt.xlsx")
