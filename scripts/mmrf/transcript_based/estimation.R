#!/usr/bin/r

## file: estimation.R
## last update: 09-09-2025


library(data.table)
library(ggbreak)
library(ggplot2)
library(readxl)
library(tidyverse)

source("C:/Users/Dell/Desktop/git_projects/TP53/scripts/fun/utils.R")

computeProb <- function(r, n){
  num <- r^abs(n - 4)
  den <- ((1 + r) ^ n) * ((1 + r) ^ (4 - n))
  #num <- (1/(1 + r))^n
  #den <- (r/(1 + r))^(4 - n)
  return(num/ den)
}

roundSig <- function(x, digits = 3) {
  return(x %/% 1 + signif(x %% 1, digits))
} 

# ---- Main ----
# setting the env
wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")
filePath <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA22/expression_estimates_transcript_based/"

# loading inputs
mmrf_counts <- fread(paste0(filePath, "MMRF_CoMMpass_IA22_salmon_transcriptUnstrandedIgFiltered_counts.tsv")) %>%
  rename_with(~"transcript", "Transcript")
mmrf_cln_per_pt <- fread(paste0("clinical_data/mmrf_cln_per_pt.txt"))

IDs <- read_excel("mmrf_tp53.xlsx", sheet = "isoforms") %>% 
  filter(NAME %in% c("p53a", "p53b", "D40p53a", "D133p53a")) %>%
  pull(`ENSEMBL mRNA`)

# filtering counts 
mmrf_counts_per_pt <- filterMMRF(mmrf_counts, mmrf_cln_per_pt) %>%
  mutate(transcript = mmrf_counts$transcript, .before = MMRF_1021)

mmrf_tp53_counts_per_pt <- as.data.frame(t(mmrf_counts_per_pt %>%
  filter(transcript %in% IDs)  %>%
  mutate(transcript = c("p53_FL", "p53β", "Δ40p53α", "Δ133p53α")) %>%
  column_to_rownames("transcript"))) 

write_tsv(mmrf_tp53_counts_per_pt, "transcript_based/mmrf_tp53_counts_per_pt.txt")


# ---- Computing ratio r and probability p: counts ----
# adding pseudocount of 1
mmrf_tp53_counts_per_pt <- mmrf_tp53_counts_per_pt + 1 #need to check if it makes sense

mmrf_tp53_prob_per_pt <- mmrf_tp53_counts_per_pt %>%
  mutate(r_Δ40_FL = roundSig(Δ40p53α/p53_FL), 
         r_Δ133_FL = roundSig(Δ133p53α/p53_FL),
         P4_FL_Δ40 = roundSig(computeProb(p53_FL, 4)),
         P3_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 3)), 
         P2_FL_Δ40 = roundSig(6 * computeProb(r_Δ40_FL, 2)), 
         P1_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 1)), 
         P0_FL_Δ40 = roundSig(computeProb(r_Δ40_FL, 0)), 
         P4_FL_Δ133 = roundSig(computeProb(p53_FL, 4)),
         P3_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 3)), 
         P2_FL_Δ133 = roundSig(6 * computeProb(r_Δ133_FL, 2)), 
         P1_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 1)), 
         P0_FL_Δ133 = roundSig(computeProb(r_Δ133_FL, 0)), 
         r_p53β_FL = roundSig(p53β/p53_FL), 
         r_p53β_Δ40 = roundSig(p53β/Δ40p53α), 
         r_p53β_Δ133 = roundSig(p53β/Δ133p53α)) %>%
  rownames_to_column("PUBLIC_ID")

write_tsv(mmrf_tp53_prob_per_pt, "transcript_based/mmrf_tp53_prob_per_pt.txt")


# ---- Computing ratio r and probability p: normalized ----
mmrf_tp53_N_per_pt <- fread("transcript_based/mmrf_tp53_per_pt.txt") %>%
  select(p53_FL, p53β, Δ40p53α, Δ133p53α)
#mmrf_tp53_N_per_pt <- mmrf_tp53_N_per_pt + 1
  
mmrf_tp53_norm_prob_per_pt <- mmrf_tp53_N_per_pt %>%
  mutate(r_Δ40_FL = roundSig(Δ40p53α/p53_FL), 
         r_Δ133_FL = roundSig(Δ133p53α/p53_FL),
         P4_FL_Δ40 = roundSig(computeProb(p53_FL, 4)),
         P3_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 3)), 
         P2_FL_Δ40 = roundSig(6 * computeProb(r_Δ40_FL, 2)), 
         P1_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 1)), 
         P0_FL_Δ40 = roundSig(computeProb(r_Δ40_FL, 0)), 
         P4_FL_Δ133 = roundSig(computeProb(p53_FL, 4)),
         P3_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 3)), 
         P2_FL_Δ133 = roundSig(6 * computeProb(r_Δ133_FL, 2)), 
         P1_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 1)), 
         P0_FL_Δ133 = roundSig(computeProb(r_Δ133_FL, 0)), 
         r_p53β_FL = roundSig(p53β/p53_FL), 
         r_p53β_Δ40 = roundSig(p53β/Δ40p53α), 
         r_p53β_Δ133 = roundSig(p53β/Δ133p53α))

#write_tsv(mmrf_tp53_prob_per_pt, "transcript_based/mmrf_tp53_prob_per_pt.txt")


# ---- Plotting: counts ----
#D40
prob_Δ40_FL <- mmrf_tp53_prob_per_pt %>% 
  select(PUBLIC_ID, P4_FL_Δ40, P3_FL_Δ40, P2_FL_Δ40, P1_FL_Δ40, P0_FL_Δ40) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
prob_Δ40_FL$group <- factor(prob_Δ40_FL$group, levels = unique(prob_Δ40_FL$group))

prob_Δ40_FL %>%
  ggplot(aes(x = group, y = probability)) +
  geom_boxplot( outlier.shape=NA )
    geom_point(aes(colour = group)) +
    stat_summary(aes(x = group, y = probability), 
                 fun.min = function(z) { quantile(z, 0.25) }, 
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.1) +
    stat_summary(aes(x = group, y = probability), 
                 fun = median, geom = "point", position = position_dodge(width = 0.5)) +
    theme_minimal()

#D133
prob_Δ133_FL <- mmrf_tp53_ratio_per_pt %>%
  select(PUBLIC_ID, P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
prob_Δ133_FL$group <- factor(prob_Δ133_FL$group, levels = unique(prob_Δ133_FL$group))

prob_Δ133_FL %>%
  ggplot(aes(x = group, y = probability)) +
    geom_point(aes(colour = group)) +
    scale_y_continuous(limits = c(0, 1)) +
    stat_summary(aes(x = group, y = probability), 
                 fun.min = function(z) { quantile(z, 0.25) }, 
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.1) +
    stat_summary(aes(x = group, y = probability), 
                 fun = median, geom = "point", position = position_dodge(width = 0.5)) +
    theme_minimal()
    
# ---- Plotting: normalized ----
#D40
norm_prob_Δ40_FL <- mmrf_tp53_norm_prob_per_pt %>% 
  select(P4_FL_Δ40, P3_FL_Δ40, P2_FL_Δ40, P1_FL_Δ40, P0_FL_Δ40) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
norm_prob_Δ40_FL$group <- factor(norm_prob_Δ40_FL$group, levels = unique(norm_prob_Δ40_FL$group))
    
norm_prob_Δ40_FL %>%
  filter(probability!=0) %>%
  ggplot(aes(x = group, y = probability, color = group)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    scale_color_manual(values = c("#2ecc71", "#3498db", "#e67e22", "#9b59b6", "#e74c3c")) + 
    scale_y_continuous(limits = c(0, 1)) +
    stat_summary(aes(x = group, y = probability, color = group), 
                 fun.min = function(z) { quantile(z, 0.25) },
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.1) +
    stat_summary(aes(x = group, y = probability, color = group), 
                 fun = median, geom = "errorbar", color = "black") +
    theme_minimal()
    
#D133
norm_prob_Δ133_FL <- mmrf_tp53_norm_prob_per_pt %>%
  select(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
norm_prob_Δ133_FL$group <- factor(norm_prob_Δ133_FL$group, levels = unique(norm_prob_Δ133_FL$group))

norm_prob_Δ133_FL %>%
  filter(!is.na(probability)) %>%
  ggplot(aes(x = group, y = probability, colour = group)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    scale_color_manual(values = c("#2ecc71", "#3498db", "#e67e22", "#9b59b6", "#e74c3c")) + 
    stat_summary(aes(x = group, y = probability), 
                 fun.min = function(z) { quantile(z, 0.25) },
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.4) +
    stat_summary(aes(x = group, y = probability), 
                 fun = median, geom = "crossbar", color = "black", width = 0.4) +
    theme_minimal()


#B
ratio_βs <- mmrf_tp53_norm_prob_per_pt %>%
  select(r_p53β_FL, r_p53β_Δ40, r_p53β_Δ133) %>%
  pivot_longer(cols = starts_with("r"), names_to = "group", values_to = "probability")
ratio_βs$group <- factor(ratio_βs$group, levels = unique(ratio_βs$group))

ratio_βs %>%
  ggplot(aes(x = group, y = probability, colour = group)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    scale_color_manual(values = c("#2ecc71", "#3498db", "#e67e22")) + 
    stat_summary(aes(x = group, y = probability), 
                 fun.min = function(z) { quantile(z, 0.25) },
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.4) +
    stat_summary(aes(x = group, y = probability), 
                 fun = median, geom = "crossbar", color = "black", width = 0.4) +
    theme_minimal()
    
### ---- WB data ----
wb_df <- read_excel("../original_db/Tabella pz TP53 (5) - label pulite.xlsx", sheet = "Densitometry") %>%
  select(p53_FL_53Kd_cont, p53_beta_gamma_45Kd_cont, delta40_p53_alpha_42Kd_cont, delta133_p53_alpha_35Kd_cont) %>%
  filter(!p53_FL_53Kd_cont=="nv") %>%
  mutate(across(everything(),  .fns = as.numeric))

wb_tp53_prob_per_pt <- wb_df %>%
  mutate(r_Δ40_FL = roundSig(delta40_p53_alpha_42Kd_cont/p53_FL_53Kd_cont), 
         r_Δ133_FL = roundSig(delta133_p53_alpha_35Kd_cont/p53_FL_53Kd_cont),
         P4_FL_Δ40 = roundSig(computeProb(p53_FL_53Kd_cont, 4)),
         P3_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 3)), 
         P2_FL_Δ40 = roundSig(6 * computeProb(r_Δ40_FL, 2)), 
         P1_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 1)), 
         P0_FL_Δ40 = roundSig(computeProb(r_Δ40_FL, 0)), 
         P4_FL_Δ133 = roundSig(computeProb(p53_FL_53Kd_cont, 4)),
         P3_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 3)), 
         P2_FL_Δ133 = roundSig(6 * computeProb(r_Δ133_FL, 2)), 
         P1_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 1)), 
         P0_FL_Δ133 = roundSig(computeProb(r_Δ133_FL, 0)), 
         r_p53β_FL = roundSig(p53_beta_gamma_45Kd_cont/p53_FL_53Kd_cont), 
         r_p53β_Δ40 = roundSig(p53_beta_gamma_45Kd_cont/delta40_p53_alpha_42Kd_cont), 
         r_p53β_Δ133 = roundSig(p53_beta_gamma_45Kd_cont/delta133_p53_alpha_35Kd_cont))

#D40
wb_prob_Δ40_FL <- wb_tp53_prob_per_pt %>% 
  select(P4_FL_Δ40, P3_FL_Δ40, P2_FL_Δ40, P1_FL_Δ40, P0_FL_Δ40) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
wb_prob_Δ40_FL$group <- factor(wb_prob_Δ40_FL$group, levels = unique(wb_prob_Δ40_FL$group))

wb_prob_Δ40_FL %>%
  ggplot(aes(x = group, y = probability, color = group)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    scale_color_manual(values = c("#2ecc71", "#3498db", "#e67e22", "#9b59b6", "#e74c3c")) + 
    stat_summary(aes(x = group, y = probability, color = group), 
                 fun.min = function(z) { quantile(z, 0.25) },
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.2) +
    stat_summary(aes(x = group, y = probability, color = group), 
                 fun = median, geom = "crossbar", color = "black", width = 0.1) +
    theme_minimal()

#D133
wb_prob_Δ133_FL <- wb_tp53_prob_per_pt %>%
  select(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
wb_prob_Δ133_FL$group <- factor(wb_prob_Δ133_FL$group, levels = unique(wb_prob_Δ133_FL$group))

wb_prob_Δ133_FL %>%
  ggplot(aes(x = group, y = probability, colour = group)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#2ecc71", "#3498db", "#e67e22", "#9b59b6", "#e74c3c")) + 
  stat_summary(aes(x = group, y = probability), 
               fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               geom = "errorbar", color = "black", size = 1.2, width = 0.2) +
  stat_summary(aes(x = group, y = probability), 
               fun = median, geom = "crossbar", color = "black", width = 0.1) +
  theme_minimal()


#B
wb_ratio_βs <- wb_tp53_prob_per_pt %>%
  select(r_p53β_FL, r_p53β_Δ40, r_p53β_Δ133) %>%
  pivot_longer(cols = starts_with("r"), names_to = "group", values_to = "probability")
wb_ratio_βs$group <- factor(wb_ratio_βs$group, levels = unique(wb_ratio_βs$group))

wb_ratio_βs %>%
  ggplot(aes(x = group, y = probability, colour = group)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    scale_color_manual(values = c("#2ecc71", "#3498db", "#e67e22")) + 
    stat_summary(aes(x = group, y = probability), 
                 fun.min = function(z) { quantile(z, 0.25) },
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.2) +
    stat_summary(aes(x = group, y = probability), 
                 fun = median, geom = "crossbar", color = "black", width = 0.1) +
    theme_minimal() +
  scale_y_break(c(2000, 4000))
