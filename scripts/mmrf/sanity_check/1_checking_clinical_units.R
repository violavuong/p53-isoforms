# !/usr/bin/r


# file: 1_checking_clinical_units
# aim: check lab-clinical data redundancy across releases IA22 and IA24 to figure out units
# last update: 06-08-2026


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
# setting env
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/")
mmrfDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/clinical_data/labs_deid.csv"
outDir <- "/mmrf/"


# ---- Cohort ----
mmrf_tp53_cohort_biochem <- read.xlsx("db_clinico_e_genomico_gennaio2026.xlsx") %>%
  select(PUBLIC_ID, ALB, CREATININE, HE, LDH, LDH_level, serum_M_prot, PLT, serum_B2_mglob) %>%
  arrange(PUBLIC_ID)


# ---- Labs IA24 ----
mmrf_IA24_labs <- fread(paste0(mmrfDir, "labs_deid.csv")) %>%
  mutate(PUBLIC_ID = toupper(public_id)) %>%
  filter(PUBLIC_ID %in% mmrf_tp53_cohort_biochem$PUBLIC_ID, lab_test %in% c("hemoglobin", "platelets", "serum_albumin", "serum_creatinine", "serum_ldh", "b2m"), baseline_lab=="yes") %>%
  select(PUBLIC_ID, lab_test, lab_test_result, lab_test_result_unit) %>%
  distinct(PUBLIC_ID, lab_test, .keep_all = TRUE) %>%
  pivot_wider(id_cols = PUBLIC_ID, names_from = lab_test, values_from = c(lab_test_result, lab_test_result_unit))
  #select(public_id, b2m=lab_test_result_b2m, b2m_unit=lab_test_result_unit_b2m, serum_creatinine=lab_test_result_serum_creatinine, serum_creatinine_unit=lab_test_result_unit_serum_creatinine)


# ---- Merge ----
mmrf_labs_IA22_IA24 <- left_join(mmrf_tp53_cohort_biochem, mmrf_IA24_labs, by = "PUBLIC_ID") %>%
  select(PUBLIC_ID, 
         IA22_He=HE, IA22_He=lab_test_result_hemoglobin, He_unit=lab_test_result_unit_hemoglobin, 
         IA22_plt=PLT, IA24_plt=lab_test_result_platelets, plt_unit=lab_test_result_unit_platelets, 
         IA22_b2m=serum_B2_mglob, IA24_b2m=lab_test_result_b2m, b2m_unit=lab_test_result_unit_b2m, 
         IA22_alb=ALB, IA24_alb=lab_test_result_serum_albumin, alb_unit=lab_test_result_unit_serum_albumin, 
         IA22_creatinine=CREATININE, IA24_creatinine=lab_test_result_serum_creatinine, creatinine_unit=lab_test_result_unit_serum_creatinine, 
         IA22_LDH=LDH, IA24_LDH=lab_test_result_serum_ldh, LDH_unit=lab_test_result_unit_serum_ldh,
         IA22_LDH_level=LDH_level)

View(mmrf_labs_IA22_IA24 %>% select(PUBLIC_ID, IA22_creatinine, IA24_creatinine, creatinine_unit))


# ---- Check He ----
mmrf_IA24_He <- fread(paste0(mmrfDir, "labs_deid.csv")) %>%
  mutate(PUBLIC_ID = toupper(public_id)) %>%
  filter(PUBLIC_ID %in% mmrf_tp53_cohort_biochem$PUBLIC_ID, lab_test=="hemoglobin")
