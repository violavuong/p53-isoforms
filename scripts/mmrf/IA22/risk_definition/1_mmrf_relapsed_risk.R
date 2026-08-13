# !/usr/bin/r


# file: 1_mmrf_relapsed_risk
# aim: defining risk groups in mmrf relapsed pts
# last update: 06-08-26


library(data.table)
library(openxlsx)
library(tidyverse)
library(waldo)


# ---- Main ----
# setting env 
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")


# ---- Cohort ----
mmrf_relapse_per_pt <- fread("mmrf_relapsed_cohort.txt") %>%
  rename_at(c(53:59), ~c("p53_FL", "p53_beta", "p53_gamma", "delta133_p53_gamma", "delta133_p53_alpha", "delta133_p53_beta", "delta40_p53_alpha"))

# ---- Exploration pt-specific ----
dup_pts <- mmrf_relapse_per_pt[duplicated(mmrf_relapse_per_pt$PUBLIC_ID) | duplicated(mmrf_relapse_per_pt$PUBLIC_ID, fromLast = TRUE), ] 

clean_dup_pts <- dup_pts %>%
  arrange(PUBLIC_ID, desc(therendy), desc(nchar(thercat))) %>%
  distinct(PUBLIC_ID, .keep_all = TRUE)

clean_dup_pts <- rbind(clean_dup_pts, dup_pts %>% filter(SPECTRUM_SEQ %in% c("MMRF_1790_2", "MMRF_1790_3")))


# ---- Clean the cohort ----
mmrf_relapsed <- rbind(mmrf_relapse_per_pt %>% filter(!PUBLIC_ID %in% dup_pts$PUBLIC_ID), clean_dup_pts) %>%
  arrange(PUBLIC_ID) %>%
  filter(!is.na(`1p`))


# ---- Alteration classification ----
median_relapsed_per_tp53 <- mmrf_relapsed %>% 
  summarise(across(c(53:59), ~median(.[. > 0])))

mmrf_relapsed_cat <- mmrf_relapsed %>% 
  mutate(across(contains("p53", ignore.case = FALSE), ~ factor(ifelse(.x==0, "absent", ifelse(.x > median_relapsed_per_tp53[[cur_column()]], "high", "low")), levels = c("absent", "low", "high")), .names = "{.col}_exp")) %>%
  mutate(across(ends_with("exp"), ~ifelse(.=="absent", 1, 0), .names = "{col}_absent")) %>%
  mutate(across(ends_with("exp"), ~ifelse(.=="low", 1, 0), .names = "{col}_low")) %>%
  mutate(across(ends_with("exp"), ~ifelse(.=="high", 1, 0), .names = "{col}_high")) %>%
  mutate(resp_group = ifelse(bestrespsh=="", "nv", ifelse(bestrespsh %in% c("sCR", "CR", "VGPR"), "responder", "non_responder"))) %>%
  mutate(p53FL_p53beta_balanced = ifelse(p53_FL_exp==p53_beta_exp, 1, 0)) %>%
  mutate(across(c(60:104), ~ifelse(. > 2.10, "amp", ifelse(. < 1.90, "del", "normal")), .names= "label_{col}_alt")) %>%
  mutate(across(c(60:104), ~ifelse(. > 2.10, 1, 0), .names= "label_{col}_amp")) %>%
  mutate(across(c(60:104), ~ifelse(. < 1.90, 1, 0), .names= "label_{col}_del")) %>%
  mutate(fish_t_4_14 = ifelse(NSD2_CALL==1 | str_detect(t_IGH, "4;14"), 1, 0), 
         fish_t_6_14 = ifelse(CCND3_CALL==1 | str_detect(t_IGH, "6;14"), 1, 0),
         fish_t_8_14 = ifelse(MYC_CALL==1 | str_detect(t_IGH, "8;14"), 1, 0),
         fish_t_8_14_2 = ifelse(MAFA_CALL==1 | str_detect(t_IGH, "8;14"), 1, 0),
         fish_t_11_14 = ifelse(CCND1_CALL==1 | str_detect(t_IGH, "11;14"), 1, 0),
         fish_t_12_14 = ifelse(CCND2_CALL==1 | str_detect(t_IGH, "12;14"), 1, 0),
         fish_t_14_16 = ifelse(MAF_CALL==1 | str_detect(t_IGH, "14;16"), 1, 0),
         fish_t_14_20 = ifelse(MAFB_CALL==1 | str_detect(t_IGH, "14;20"), 1, 0)) %>%
  mutate(mutTP53_POS = ifelse(!is.na(TP53_POS), 1, 0)) %>%
  mutate(label_3_amp = ifelse(label_3p_amp==1 | label_3q_amp==1, 1, 0),
         label_5_amp = ifelse(label_5p_amp==1 | label_5q_amp==1, 1, 0),
         label_7_amp = ifelse(label_7p_amp==1 | label_7q_amp==1, 1, 0),
         label_9_amp = ifelse(label_9p_amp==1 | label_9q_amp==1, 1, 0),
         label_11_amp = ifelse(label_11p_amp==1 | label_11q_amp==1, 1, 0),
         label_19_amp = ifelse(label_19p_amp==1 | label_19q_amp==1, 1, 0),
         label_21_amp = ifelse(label_21p_amp==1 | label_21q_amp==1, 1, 0), 
         tot_amp = label_3_amp + label_5_amp + label_7_amp + label_9_amp + label_11_amp + label_15q_amp + label_19q_amp, label_21q_amp,
         hd = ifelse(tot_amp > 2, 1, 0), 
         hd_with_19q = ifelse(hd==1 & label_19q_amp==1, 1, 0), 
         hd_with_21q = ifelse(hd==1 & label_21q_amp==1, 1, 0))
  

write.xlsx(mmrf_relapsed_cat, "../db_clinico_genomico_relapse_300726.xlsx")


# ---- Risk Group ----
mmrf_relapsed_cat <- read.xlsx("../db_clinico_genomico_relapse_300726.xlsx")

mmrf_relapsed_risk <- mmrf_relapsed_cat %>%
  mutate(risk_group = ifelse((p53_FL_exp_high==1 & p53_beta_exp_absent==0 & delta133_p53_beta_exp_absent==1) & (hd_with_21q==1 | label_11p_amp==1 | label_18p_amp==1 | label_18q_amp==1), "favorable",
         ifelse(((p53_FL_exp_high==0 & p53_beta_exp_high==0) | p53FL_p53beta_balanced==0 | delta133_p53_beta_exp_absent==0 | p53_gamma_exp_high==0) & (label_11p_del==1 | fish_t_8_14_2==1), "poor", "other"))) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, p53_FL_exp_high, p53_beta_exp_absent, delta133_p53_beta_exp_absent, hd_with_21q, label_11p_amp, label_18p_amp, label_18q_amp, 
         p53_beta_exp_high, p53FL_p53beta_balanced, delta133_p53_beta_exp_absent, p53_gamma_exp_high, label_11p_del, fish_t_8_14_2, risk_group)


table(mmrf_relapsed_risk$risk_group)

write.xlsx(mmrf_relapsed_risk, "../mmrf_relapsed_risk_groups_300726.xlsx")


# ---- Adding baseline group risk ----
# genomic/transcriptomic profile at baseline
mmrf_baseline_per_pt <- read.xlsx("../db_clinico_genomico_sopravvivenza_1609257.xlsx") %>%
  filter(PUBLIC_ID %in% mmrf_relapsed_risk$PUBLIC_ID) %>%
  mutate(SPECTRUM_SEQ = paste0(PUBLIC_ID, "_1"), 
         line = 1, 
         TP53_MUT = ifelse(is.na(TP53_POS), 0, 1)) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, line, 
         resp_sh, resp_group, 
         start_dy=start_line_1, end_dy=end_line_1, best_respdy=best_resp_dy_line_1, 
         PFS_dy=PFS_days, PFS_event, PFS_time=PFS_date, 
         OS_date, OS_event, OS_time=OS_days, 
         p53_FL=p53_FL_cont, p53_FL_exp=p53_FL, 
         p53_beta=p53_beta_cont, p53_beta_exp=p53_beta, 
         p53_gamma=p53_gamma_cont, p53_gamma_exp=p53_gamma, 
         delta40_p53_alpha=delta40_p53_alpha_cont, delta40_p53_alpha_exp=delta40_p53_alpha, 
         delta133_p53_alpha=delta133_p53_alpha_cont, delta133_p53_alpha_exp=delta133_p53_alpha, 
         delta133_p53_beta=delta133_p53_gamma_cont, delta133_p53_beta_exp=delta133_p53_beta, 
         delta133_p53_gamma=delta133_p53_gamma_cont, delta133_p53_gamma_exp=delta133_p53_gamma, 
         66:167, TP53_MUT)

mmrf_baseline_risk <- read.xlsx("../risk_groups/mmrf_tp53_risk_groups_210726.xlsx") %>%
  filter(PUBLIC_ID %in% mmrf_relapsed_risk$PUBLIC_ID) 
  
mmrf_E_per_pt <- left_join(mmrf_baseline_per_pt, mmrf_baseline_risk, by = "PUBLIC_ID") %>%
  arrange(PUBLIC_ID)

#write.xlsx(mmrf_E_per_pt, "../relapse/mmrf_E_per_relapsed_pt_060826.xlsx")  # manual revision in excel (changed some cols name), reloading 
mmrf_E_per_pt <- read.xlsx("../relapse/mmrf_E_per_relapsed_pt_060826.xlsx")

# genomic/transcriptomic profile at relapse
mmrf_R_per_pt <- mmrf_relapsed_cat %>%
  mutate(TP53_MUT = ifelse(is.na(TP53_POS), 0, 1)) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, line, resp_sh=bestrespsh, resp_group, 
         start_dy, end_dy, best_respdy, 
         PFS_dy, PFS_event, PFS_time, 
         OS_date, OS_event, OS_time,
         p53_FL, p53_FL_exp, 
         p53_beta, p53_beta_exp, 
         p53_gamma, p53_gamma_exp,
         delta40_p53_alpha, delta40_p53_alpha_exp, 
         delta133_p53_alpha, delta133_p53_alpha_exp, 
         delta133_p53_beta, delta133_p53_beta_exp, 
         delta133_p53_gamma, delta133_p53_gamma_exp,
         60:113, TP53_MUT, 
         starts_with("label")) %>%
  left_join(mmrf_relapsed_risk, by = c("PUBLIC_ID", "SPECTRUM_SEQ"))

#write.xlsx(mmrf_R_per_pt, "../relapse/mmrf_R_per_relapsed_pt_060826.xlsx")  # manual revision in excel (changed/removed some cols name), reloading 
mmrf_R_per_pt <- read.xlsx("../relapse/mmrf_R_per_relapsed_pt_060826.xlsx")
colnames(mmrf_R_per_pt)==colnames(mmrf_E_per_pt)


# merge
mmrf_E_R_per_pt <- rbind(mmrf_E_per_pt, mmrf_R_per_pt) %>%
  arrange(SPECTRUM_SEQ)

write.xlsx(mmrf_E_R_per_pt, "../relapse/mmrf_E_R_per_relapsed_pt_060826.xlsx")
