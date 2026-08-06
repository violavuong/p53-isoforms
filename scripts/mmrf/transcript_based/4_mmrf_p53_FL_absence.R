#!/usr/bin/r


# file: 4_mmrf_p53_FL_absence
# aim: checking if pts with p53 FL absent are also del17p
# last update: 06-08-2026


library(data.table)
library(ggplot2)
library(gtools)
library(scales)
library(tidyverse)


# ---- Main ----
# environment setting 
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")

mmrf_tp53_per_pt <- fread(paste0(wd, "/transcript_based/mmrf_tp53_per_pt.txt")) #transcriptomic
mmrf_genomic_per_pt <- fread(paste0(wd, "/genomic_based/mmrf_genomic_per_pt_class_10.txt")) #genomic


# ---- Defining two classes: FL absent / FL present ----
## TP53FL absent class 
mmrf_absent_FL_per_pt <- mmrf_tp53_per_pt %>% filter(p53_FL_exp=="absent")
mmrf_absent_profile_per_pt <- left_join(mmrf_absent_FL_per_pt, 
                                        mmrf_genomic_per_pt %>% filter(PUBLIC_ID %in% mmrf_absent_FL_per_pt$PUBLIC_ID), 
                                        by = "PUBLIC_ID")

## TP53FL present class 
mmrf_present_per_pt <- mmrf_tp53_per_pt %>% filter(p53_FL_exp!="absent")
mmrf_present_profile_per_pt <-  left_join(mmrf_present_per_pt, 
                                          mmrf_genomic_per_pt %>% filter(PUBLIC_ID %in% mmrf_present_per_pt$PUBLIC_ID), 
                                          by = "PUBLIC_ID")


# ---- Exploration ----
## frequency tables
absent_freq_tbl <- mmrf_absent_profile_per_pt %>%
  select(ends_with("alt", ignore.case = FALSE)) %>%
  pivot_longer(cols = everything(), names_to = "chrarm", values_to = "alt") %>%
  count(chrarm, alt, name = "count") %>%
  mutate(proportion = count / 32, 
         group = "FL_exp_absent") 

present_freq_tbl <- mmrf_present_profile_per_pt %>%
  select(ends_with("alt", ignore.case = FALSE)) %>%
  pivot_longer(cols = everything(), names_to = "chrarm", values_to = "alt") %>%
  count(chrarm, alt, name = "count") %>%
  mutate(proportion = count / 627, 
         group = "FL_exp_present") 


freq_tbl <- rbind(absent_freq_tbl, present_freq_tbl)
freq_tbl$chrarm <- factor(freq_tbl$chrarm, levels = mixedsort(unique(freq_tbl$chrarm)))

# plotting
freq_tbl %>%
  ggplot(aes(x = chrarm, y = proportion, fill = alt)) + 
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~group, ncol = 1) +
    scale_fill_manual(values = c("amp" = "#84B5FD", "del" = "lightskyblue", "normal" = "navajowhite")) +
    scale_y_continuous(labels = percent) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom")
    

