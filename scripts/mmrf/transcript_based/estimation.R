#!/usr/bin/r

## file: estimation.R
## last update: 23-07-2025


library(data.table)
library(ggplot2)
library(readxl)
library(tidyverse)

source("C:/Users/Dell/Desktop/git_projects/TP53/scripts/fun/utils.R")

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


# ---- Computing ratio r and probability p ----
# adding pseudocount of 1
mmrf_tp53_counts_per_pt <- mmrf_tp53_counts_per_pt + 1 #need to check if it makes sense

computeProb <- function(r, n){
  num <- r^abs(n - 4)
  den <- ((1 + r) ^ n) * ((1 + r) ^ (4 - n))
  return(num/den)
}

roundSig <- function(x, digits = 3) {
  return(x %/% 1 + signif(x %% 1, digits))
} 

mmrf_tp53_ratio_per_pt <- mmrf_tp53_counts_per_pt %>%
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

write_tsv(mmrf_tp53_ratio_per_pt, "transcript_based/mmrf_tp53_ratio_per_pt.txt")


# ---- Plotting ----
#D40
prob_Δ40_FL <- mmrf_tp53_ratio_per_pt %>% select(PUBLIC_ID, P4_FL_Δ40, P3_FL_Δ40, P2_FL_Δ40, P1_FL_Δ40, P0_FL_Δ40)
stats <- unlist(lapply(seq_along(prob_Δ40_FL)[-1], function(x) quantile(prob_Δ40_FL[, x], probs = c(0.25, 0.50, 0.75))))
Δ40_FL_stats <- data.frame(group = colnames(prob_Δ40_FL)[-1],
                           IQR1 = stats[c(1,4,7,10,13)],
                           median = stats[c(2,5,8,11,14)], 
                           IQR3 = stats[c(3,6,9,12,15)])

prob_Δ40_FL_long <- pivot_longer(prob_Δ40_FL, cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability") %>%
  left_join(Δ40_FL_stats, by = "group")
prob_Δ40_FL_long$group <- factor(prob_Δ40_FL_long$group, levels = unique(prob_Δ40_FL_long$group))

prob_Δ40_FL_long %>%
  ggplot(aes(x = group, y = probability)) +
    geom_point(aes(colour = group)) +
    #geom_errorbar(aes(ymin = IQR1, ymax = IQR3)) +
    stat_summary(fun.data = prob_Δ40_FL_long, fun.args = list(IQR1, IQR3), geom = "errorbar", color = "black", width = 0.4) +
    stat_summary(fun.y = median, geom = "point", color = "black") +
    theme_minimal()

#D133
prob_Δ133_FL <- mmrf_tp53_ratio_per_pt %>%
  select(PUBLIC_ID, P4_FL_Δ133, P3_FL_Δ133, P2_FL_Δ133, P1_FL_Δ133, P0_FL_Δ133) %>%
  pivot_longer(cols = starts_with(c("P4", "P3", "P2", "P1", "P0")), names_to = "group", values_to = "probability")
prob_Δ133_FL$group <- factor(prob_Δ133_FL$group, levels = unique(prob_Δ133_FL$group))

prob_Δ133_FL %>%
  ggplot(aes(x = group, y = probability)) +
  geom_boxplot()

