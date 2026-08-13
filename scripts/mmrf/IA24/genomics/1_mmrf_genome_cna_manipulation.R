# !/usr/bin/r


# file: 1_mmrf_genome_cna_manipulation
# aim: manipulate genome cna data to obtain a dataset per pt, per line, with a unique CN per chrarm
# last update: 12-08-26


library(data.table)
library(GenomicRanges)
library(gtools)
library(openxlsx)
library(plyranges)
library(tidyverse)


# ---- Main ----
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/genomics/")
mmrfDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/copy_number/"

libDir <- "C:/Users/violameixian.vuong2/Desktop/git-projects/unique-molecular-assay/UMA_lib/"
refGen <- "hg38"

# ---- Loading custom function ----
source("C:/Users/violameixian.vuong2/Desktop/git-projects/unique-molecular-assay/UMA_lib/fun/Popeye2.R")


# ---- Dataset Manipulation ----
mmrf_gatk_genome_seg <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA24_genome_gatk_cna.seg")) %>% 
  mutate(public_id = tolower(str_extract(SAMPLE, "MMRF_[0-9]+")), 
         line = str_extract(SAMPLE, "(?<=MMRF_[0-9]{4}_)[0-9]+"), 
         specimen = str_extract(SAMPLE, "(BM|PB)_CD138pos") %>% str_remove("_CD138pos"), 
         chr = str_remove(Chromosome, "chr"), 
         CN = 2^(Segment_Mean + 1)) %>%
  select(ID = public_id, line, specimen, chr, start = Start, end = End, CN)
 
# run Popeye2
mmrf_gatk_genome_chrarm <- Popeye2(mmrf_gatk_genome_seg, refGen, removeXY = FALSE)


# weighted mean of segments
mmrf_gatk_genome_per_chrarm <- mmrf_gatk_genome_chrarm %>%
  group_by(ID, line, specimen, chrarm) %>%
  summarise(weighted_mean_CN = (sum(CN * width)) / sum(width), 
            start = min(start), 
            end = max(end), 
            width = end - start) 

mmrf_gatk_genome_per_chrarm$chrarm <- factor(mmrf_gatk_genome_per_chrarm$chrarm, levels = unique(mixedsort(mmrf_gatk_genome_per_chrarm$chrarm)))

# ---- Transpose data ----
mmrf_gatk_genome_per_pt <- mmrf_gatk_genome_per_chrarm %>%
  dplyr::rename(public_id = ID) %>%
  pivot_wider(id_cols = c(public_id, line, specimen), 
              names_from = chrarm, 
              values_from = weighted_mean_CN, 
              names_sort = TRUE) %>%
  arrange(public_id, line)

names(mmrf_gatk_genome_per_pt)[4:48] <- paste0("chrarm_", names(mmrf_gatk_genome_per_pt)[4:48])

write.xlsx(mmrf_gatk_genome_per_pt, "/complete_db/mmrf_gatk_genome_per_pt_ia24.xlsx")
