# !/usr/bin/r


# file: 3_mmrf_risk_calculation
# aim: compute the risk groups of first line pts (n=659)
# last update: 07-08-2026


library(data.table)
library(openxlsx)
library(tidyverse)

library(survminer)
library(survival)
library(survivalROC)

# ---- Main ----
# setting environment
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/sanity_check/")


# ---- Cohort ----
# complete dataset with genomic, transcriptomic, clinical data, each categorized accordingly
mmrf_first_line_per_pt <- read.xlsx("mmrf_first_line_per_pt.xlsx")


# ---- Compute the risk ----
mmrf_first_with_risk_per_pt <- mmrf_first_line_per_pt %>%
  mutate(risk_group = ifelse((p53_FL_exp_high==1 & p53_beta_exp_absent==0 & delta133_p53_beta_exp_absent==1) & (hd_with_21q==1 | label_11p_alt_amp==1 | label_18p_alt_amp==1 | label_18q_alt_amp==1), "favorable",
                             ifelse(((p53_FL_exp_high==0 & p53_beta_exp_high==0) | p53FL_p53beta_balanced==0 | delta133_p53_beta_exp_absent==0 | p53_gamma_exp_high==0) & (label_11p_alt_del==1 | fish_t_8_14_2==1), 
                                    "poor", "other")))

# favorable=175 pts, other=470pts, poor=14
table(mmrf_first_with_risk_per_pt$risk_group)


# ---- Double-checking ----
# n=313 pts, median p53_FL: 12.49375
p53_FL_high_pts <- mmrf_first_with_risk_per_pt %>%
  filter(p53_FL_exp_high==1) %>%
  select(PUBLIC_ID, p53_FL, p53_FL_exp, p53_FL_exp_high)


# n=658 pts (only 1 pt does not have beta expressed)
p53_beta_expressed_pts <- mmrf_first_with_risk_per_pt %>%
  filter(p53_beta_exp_absent==0) %>%
  select(PUBLIC_ID, p53_beta, p53_beta_exp, p53_beta_exp_absent)


# n=640 pts
delta133_p53_beta_absent_pts <- mmrf_first_with_risk_per_pt %>%
  filter(delta133_p53_beta_exp_absent==1) %>%
  select(PUBLIC_ID, delta133_p53_beta, delta133_p53_beta_exp, delta133_p53_beta_exp_absent)


# n=170 pts
hd_with_21q_pts <- mmrf_first_with_risk_per_pt %>%
  filter(hd_with_21q==1) %>%
  select(PUBLIC_ID, tot_amp, label_21_amp, label_21q_alt_amp, hd, hd_with_21q)


# n=262 pts
amp_11p_pts <- mmrf_first_with_risk_per_pt %>%
  filter(label_11p_alt_amp==1) %>%
  select(PUBLIC_ID, `11p`, `11p_alt`, label_11p_alt_amp)


# n=120 pts
amp_18p_pts <- mmrf_first_with_risk_per_pt %>%
  filter(label_18p_alt_amp==1) %>%
  select(PUBLIC_ID, `18p`, `18p_alt`, label_18p_alt_amp)


# n=132 pts
amp_18q_pts <- mmrf_first_with_risk_per_pt %>%
  filter(label_18q_alt_amp==1) %>%
  select(PUBLIC_ID, `18q`, `18q_alt`, label_18q_alt_amp)


# n=346 pts (346+313=659 pts), median p53_FL: 12.49375
p53_FL_non_high_pts <- mmrf_first_with_risk_per_pt %>%
  filter(p53_FL_exp_high==0) %>%
  select(PUBLIC_ID, p53_FL, p53_FL_exp, p53_FL_exp_high)


# n=330 pts, median p53_beta: 5.191551
p53_beta_non_high_pts <- mmrf_first_with_risk_per_pt %>%
  filter(p53_beta_exp_high==0) %>%
  select(PUBLIC_ID, p53_beta, p53_beta_exp, p53_beta_exp_high)


# n=314 pts
p53FL_p53beta_non_balanced_pts <- mmrf_first_with_risk_per_pt %>%
  filter(p53FL_p53beta_balanced==0) %>%
  select(PUBLIC_ID, p53_FL_exp, p53_beta_exp, p53FL_p53beta_balanced)


# n=19 pts
delta133_p53_beta_expressed_pts <- mmrf_first_with_risk_per_pt %>%
  filter(delta133_p53_beta_exp_absent==0) %>%
  select(PUBLIC_ID, delta133_p53_beta, delta133_p53_beta_exp, delta133_p53_beta_exp_absent)


# n=564 pts
p53_gamma_non_high_pts <-  mmrf_first_with_risk_per_pt %>%
  filter(p53_gamma_exp_high==0) %>%
  select(PUBLIC_ID, p53_gamma, p53_gamma_exp, p53_gamma_exp_high)


# n=13 pts
del_11p_pts <-  mmrf_first_with_risk_per_pt %>%
  filter(label_11p_alt_del==1) %>%
  select(PUBLIC_ID, `11p`, `11p_alt`, label_11p_alt_del)


# n=4 pts
t_8_14_MAFA_pts <- mmrf_first_with_risk_per_pt %>%
  filter(fish_t_8_14_2==1) %>%
  select(PUBLIC_ID, fish_t_8_14_2)



# ---- Triple-checking ----
mmrf_risk_groups_per_pt <- mmrf_first_with_risk_per_pt %>%
  select(PUBLIC_ID, p53_FL_exp_high, p53_beta_exp_absent, p53_beta_exp_high, p53_gamma_exp_high, delta133_p53_beta_exp_absent, p53FL_p53beta_balanced, label_11p_alt_amp, label_11p_alt_del, 
         label_18p_alt_amp, label_18q_alt_amp, hd_with_21q, fish_t_8_14_2, risk_group)
  
  
