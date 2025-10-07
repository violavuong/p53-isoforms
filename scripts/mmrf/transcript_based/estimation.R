#!/usr/bin/r

## file: estimation.R
## last update: 07-10-2025


library(data.table)
library(gghalves)
library(ggplot2)
library(readxl)
library(scales)
library(tidyverse)


computeProb <- function(r, n){
  num <- r^abs(n - 4)
  den <- ((1 + r) ^ n) * ((1 + r) ^ (4 - n))
  return(num/ den)
}

roundSig <- function(x, digits = 3) {
  return(x %/% 1 + signif(x %% 1, digits))
} 

# ---- Main ----
# setting the env
wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")

# loading input
mmrf_tp53_N_per_pt <- fread("transcript_based/mmrf_tp53_per_pt.txt") %>%
  select(PUBLIC_ID, p53_FL, p53_FL_exp, Δ40p53α, Δ40p53α_exp, Δ133p53α, Δ133p53α_exp)

# ---- Computing ratio r and probability p ----
mmrf_tp53_prob_per_pt <- mmrf_tp53_N_per_pt %>%
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
         P0_FL_Δ133 = roundSig(computeProb(r_Δ133_FL, 0)))

write_tsv(mmrf_tp53_prob_per_pt, "transcript_based/mmrf_tp53_prob_per_pt.txt")


# ---- Plotting: ggplot2 ----
#D40
prob_Δ40_FL <- mmrf_tp53_prob_per_pt %>% 
  select(P4_FL_Δ40, P3_FL_Δ40, P2_FL_Δ40, P1_FL_Δ40, P0_FL_Δ40) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
prob_Δ40_FL$group <- factor(prob_Δ40_FL$group, levels = unique(prob_Δ40_FL$group))
    
prob_Δ40_FL_boxplot <- prob_Δ40_FL %>%
  filter(!is.na(probability)) %>%
  ggplot(aes(x = group, y = probability, color = group)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    labs(x = "") +
    scale_color_manual(values = c("navajowhite", "burlywood3", "lightskyblue", "#84B5FD", "navy")) + 
    stat_summary(aes(x = group, y = probability, color = group), 
                 fun.min = function(z) { quantile(z, 0.25) },
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.3) +
    stat_summary(aes(x = group, y = probability, color = group), 
                 fun = median, geom = "errorbar", color = "black", size = 0.5, width = 0.2) +
    theme_minimal()

ggsave("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/visualization/mmrf_prob_Δ40_FL.png", 
       prob_Δ40_FL_boxplot, bg = "white", dpi = 400, height = 10, width = 10)
    
#D133
prob_Δ133_FL <- mmrf_tp53_prob_per_pt %>%
  select(P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
prob_Δ133_FL$group <- factor(prob_Δ133_FL$group, levels = unique(prob_Δ133_FL$group))

prob_Δ133_FL_boxplot <- prob_Δ133_FL %>%
  filter(!is.na(probability)) %>%
  ggplot(aes(x = group, y = probability, colour = group)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    labs(x = "") +
    scale_color_manual(values = c("navajowhite", "burlywood3", "lightskyblue", "#84B5FD", "navy")) + 
    stat_summary(aes(x = group, y = probability), 
                 fun.min = function(z) { quantile(z, 0.25) },
                 fun.max = function(z) { quantile(z, 0.75) },
                 geom = "errorbar", color = "black", size = 1.2, width = 0.3) +
    stat_summary(aes(x = group, y = probability), 
                 fun = median, geom = "crossbar", color = "black", size = 0.5, width = 0.2) +
    theme_minimal()

ggsave("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/visualization/mmrf_prob_Δ133_FL.png", 
       prob_Δ133_FL_boxplot, bg = "white", dpi = 400, height = 10, width = 10)


# ---- Plotting: gghalves ----
#D40
prob_log_Δ40_FL_boxplot <- prob_Δ40_FL %>%
  filter(probability!=0) %>%
  ggplot(aes(x = group, y = probability, colour = group)) +
    geom_half_boxplot(center = TRUE, nudge = 0.01, width = 0.3) +
    geom_half_violin(aes(fill = group), side = "r") + 
    geom_half_dotplot(binaxis = "y", dotsize = 0.1, drop = TRUE, method = "histodot") +
    labs(x = "", y = "log10(probability)") +
    scale_color_manual(values = c("navajowhite", "burlywood3", "lightskyblue", "#84B5FD", "navy")) +
    scale_fill_manual(values = c("navajowhite", "burlywood3", "lightskyblue", "#84B5FD", "navy")) +
    scale_y_continuous(trans = "log10", expand = expansion(add = 1)) + #tried trans = pseudo_log_trans(sigma = 10^(-2), base = 10)
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/visualization/mmrf_norm_prob_log_Δ40_FL.png", 
       prob_log_Δ40_FL_boxplot, bg = "white", dpi = 400, height = 10, width = 20)

#D133
prob_log_Δ133_FL_boxplot <- prob_Δ133_FL %>%
  filter(probability!=0) %>%
  ggplot(aes(x = group, y = probability, colour = group)) +
    geom_half_boxplot(center = TRUE, nudge = 0.01, width = 0.3) +
    geom_half_violin(aes(fill = group), side = "r") + 
    geom_half_dotplot(binaxis = "y", binwidth = 0.25, dotsize = 0.1, drop = TRUE, method = "histodot") +
    labs(x = "", y = "log10(probability)") +
    scale_color_manual(values = c("navajowhite", "burlywood3", "lightskyblue", "#84B5FD", "navy")) +
    scale_fill_manual(values = c("navajowhite", "burlywood3", "lightskyblue", "#84B5FD", "navy")) +
    scale_y_continuous(trans = "log10", expand = expansion(add = 1)) + #tried trans = pseudo_log_trans(sigma = 10^(-2), base = 10)
    theme_minimal() +
    theme(legend.position = "bottom")


ggsave("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/visualization/mmrf_norm_prob_log_Δ133_FL.png", 
       prob_log_Δ133_FL_boxplot, bg = "white", dpi = 400, height = 10, width = 20)

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

