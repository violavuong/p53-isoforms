#!/usr/bin/r

## file: estimation.R
## last update: 22-07-2025


library(data.table)
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

mmrf_tp53_tetra <- read_excel("mmrf_tp53.xlsx", sheet = "isoforms") %>% 
  filter(NAME %in% c("p53a", "D40p53a", "D133p53a"))
IDs <- mmrf_tp53_tetra$`ENSEMBL mRNA`

# filtering counts 
mmrf_counts_per_pt <- filterMMRF(mmrf_counts, mmrf_cln_per_pt) %>%
  mutate(transcript = mmrf_counts$transcript, .before = MMRF_1021)

mmrf_tp53_counts_per_pt <- as.data.frame(t(mmrf_counts_per_pt %>%
  filter(transcript %in% IDs)  %>%
  mutate(transcript = c("p53_FL", "Δ40p53α", "Δ133p53α")) %>%
  column_to_rownames("transcript"))) 

write_tsv(mmrf_tp53_counts_per_pt, "transcript_based/mmrf_tp53_counts_per_pt.txt")


# ---- Computing ratio r and probability p ----
# adding pseudocount of 1
mmrf_tp53_counts_per_pt <- mmrf_tp53_counts_per_pt + 1 #need to check if it makes sense

computeProb <- function(r, n){
  return(r^abs(n-4) / ((1 + r) ^ n * (1 + r) ^ (4 - n)))
}

roundSig <- function(x, digits = 3) {
  return(x %/% 1 + signif(x %% 1, digits))
} 

mmrf_tp53_ratio_per_pt <- mmrf_tp53_counts_per_pt %>%
  mutate(r_Δ40_FL = roundSig(Δ40p53α/p53_FL), 
         r_Δ133_FL = roundSig(Δ133p53α/p53_FL),
         P3_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 3)), 
         P2_FL_Δ40 = roundSig(6 * computeProb(r_Δ40_FL, 2)), 
         P1_FL_Δ40 = roundSig(4 * computeProb(r_Δ40_FL, 1)), 
         P0_FL_Δ40 = roundSig(computeProb(r_Δ40_FL, 0)), 
         P3_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 3)), 
         P2_FL_Δ133 = roundSig(6 * computeProb(r_Δ133_FL, 2)), 
         P1_FL_Δ133 = roundSig(4 * computeProb(r_Δ133_FL, 1)), 
         P0_FL_Δ133 = roundSig(computeProb(r_Δ133_FL, 0))) %>%
  rownames_to_column("PUBLIC_ID")

write_tsv(mmrf_tp53_ratio_per_pt, "transcript_based/mmrf_tp53_ratio_per_pt.txt")

