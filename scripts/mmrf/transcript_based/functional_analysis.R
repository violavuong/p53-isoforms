#!/usr/bin/r

## file: functional_analysis.R
## last update: 07-10-2025


library(crosstable)
library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
# environment setting 
wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")
mmrf_tp53_prob_per_pt <- fread(paste0(wd, "/transcript_based/mmrf_tp53_prob_per_pt.txt")) #transcriptomic
mmrf_genomic_per_pt <- fread(paste0(wd, "/genomic_based/mmrf_genomic_per_pt_class_10.txt")) #genomic


# ---- Δ133 ----
# subsetting
mmrf_Δ133_prob_per_pt <- mmrf_tp53_prob_per_pt %>%
  select(PUBLIC_ID, p53_FL, p53_FL_exp, Δ133p53α, Δ133p53α_exp, ends_with("Δ133")) %>%
  mutate(tetramer = case_when(P4_FL_Δ133==pmax(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) ~ "P4_FL_Δ133", 
                              P3_FL_Δ133==pmax(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) ~ "P3_FL_Δ133", 
                              P2_FL_Δ133==pmax(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) ~ "P2_FL_Δ133", 
                              P1_FL_Δ133==pmax(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) ~ "P1_FL_Δ133", 
                              P0_FL_Δ133==pmax(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) ~ "P0_FL_Δ133"), 
        tetramer_prob = pmax(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133))


# cohort description
## how many pts fall within each tetramer category?
mmrf_Δ133_count <- mmrf_Δ133_prob_per_pt %>%
  group_by(tetramer) %>%
  summarise(count = n())

## how pts separete themselves based on tetramer, FL and 133 exp
mmrf_Δ133_count_per_class <- mmrf_Δ133_prob_per_pt %>%
  select(p53_FL_exp, Δ133p53α_exp, tetramer) 
mmrf_Δ133_count_per_class$tetramer <- factor(mmrf_Δ133_count_per_class$tetramer, levels = addNA(c("P4_FL_Δ133", "P3_FL_Δ133", "P2_FL_Δ133", "P1_FL_Δ133", "P0_FL_Δ133")))

mmrf_Δ133_tbl <- crosstable(mmrf_Δ133_count_per_class, by = c(p53_FL_exp, Δ133p53α_exp), label = FALSE, total = TRUE)
write.xlsx(mmrf_Δ133_tbl, "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/mmrf_Δ133_prob_per_class.xlsx")


# ---- Δ133: exploration by class ----
# high p53FL - low Δ133: 142 pts
high_FL_low_Δ133_per_pt <- left_join(mmrf_Δ133_prob_per_pt %>% filter(p53_FL_exp=="high", Δ133p53α_exp=="low"), 
                                     mmrf_genomic_per_pt %>% filter(PUBLIC_ID %in% high_FL_low_Δ133_per_pt$PUBLIC_ID), by = "PUBLIC_ID")

tmp <- high_FL_low_Δ133_per_pt %>%
  select(ends_with("alt", ignore.case = FALSE)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "chrarm",
    values_to = "alt"
  ) %>%
  filter(alt %in% c("normal", "amp", "del")) %>% # Keep only relevant values
  count(chrarm, alt, name = "Count")

tmp$chrarm <- factor(tmp$chrarm, levels = gtools::mixedsort(unique(tmp$chrarm)))

View(tmp %>% filter(alt!="normal"))


crosstable(tmp, label = FALSE, total = TRUE) %>% flextable::flextable()


# low p53FL - low Δ133: 119 pts
high_FL_low_Δ133_per_pt <- left_join(mmrf_Δ133_prob_per_pt %>% filter(p53_FL_exp=="high", Δ133p53α_exp=="low"), 
                                     mmrf_genomic_per_pt %>% filter(PUBLIC_ID %in% high_FL_low_Δ133_per_pt$PUBLIC_ID), by = "PUBLIC_ID")

tmp <- high_FL_low_Δ133_per_pt %>%
  select(ends_with("alt", ignore.case = FALSE)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "chrarm",
    values_to = "alt"
  ) %>%
  filter(alt %in% c("normal", "amp", "del")) %>% # Keep only relevant values
  count(chrarm, alt, name = "Count")

tmp$chrarm <- factor(tmp$chrarm, levels = gtools::mixedsort(unique(tmp$chrarm)))

View(tmp %>% filter(alt!="normal"))


crosstable(tmp, label = FALSE, total = TRUE) %>% flextable::flextable()





















