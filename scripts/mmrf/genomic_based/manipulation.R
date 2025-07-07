# !/usr/bin/r

# file: manipulation.R
# description: data extraction and handling of compelling CoMMpass genomic data
# last update: 07-07-2025

library(data.table)
library(GenomicRanges)
library(plyranges)
library(tidyverse)

source("C:/Users/Dell/Desktop/git_projects/UMA/UMA_lib/fun/Popeye2.R")

# ---- Main ----
# working directories
wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/genomic_based/")
filePath <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA22/"


# ---- GATK genome copy number ----
gene_CN_df <- fread(paste0(filePath, "copy_number/MMRF_CoMMpass_IA22_genome_gatk_cna.seg")) %>%
  filter(str_detect(SAMPLE, "1_BM_CD138pos")) %>%
  mutate(CN = 2^(Segment_Mean + 1), 
         Chromosome = str_remove(Chromosome, "chr")) %>%
  rename_with(~"ID", SAMPLE)

chrarm_CN_df <- Popeye2(gene_CN_df, "hg38", removeXY = FALSE)
chrarm_CN_df$ID <- str_remove(chrarm_CN_df$ID, "_1_BM_CD138pos")

write_tsv(chrarm_CN_df, "WGS_CN_per_chrarm.tsv")


# ---- FISH ----
fish_df <- fread(paste0(filePath, "seqFISH/MMRF_CoMMpass_IA22_genome_tumor_only_mm_igtx_pairoscope.tsv")) %>%
  filter(str_detect(SAMPLE, "1_BM_CD138pos")) %>%
  mutate(ID = str_remove(SAMPLE, "_1_BM_CD138pos")) %>%
  select(ID, ends_with("CALL"))

write_tsv(fish_df, "canonical_t_IgH_FISH.tsv")


# ---- Non-synonymous SNVs ----
ns_snv_df <- fread(paste0(filePath, "somatic_mutations/MMRF_CoMMpass_IA22_combined_vcfmerger2_All_Canonical_NS_Variants.tsv")) %>%
  filter(str_detect(SAMPLE, "1_BM_CD138pos"), 
         `ANN[*].GENE`=="TP53") %>%
  mutate(ID = str_remove(SAMPLE, "_1_BM_CD138pos"), 
         VAF = `GEN[1].AD[1]`/(`GEN[1].AD[1]` + `GEN[1].AD[0]`)) %>%
  select(ID, CHROM, POS, ID, REF, ALT, MUTECT2, STRELKA2, VARDICT, OCTOPUS, LANCET, 
         EFFECT = `ANN[*].EFFECT`, HUGO_ID = `ANN[*].GENE`, ENSG_ID = `ANN[*].GENEID`, CDS_ANN = `ANN[*].HGVS_C`, PROT_ANN = `ANN[*].HGVS_P`, 
         ALT_DP = `GEN[1].AD[1]`, REF_DP = `GEN[1].AD[0]`, VAF)

write_tsv(ns_snv_df, "tp53_SNVs.tsv")


# ---- NGS translocations ----
chr_partners <- c("chr4", "chr6", "chr11", "chr16", "chr20")

# DELLY
delly_df <- fread(paste0(filePath, "structural_event/MMRF_CoMMpass_IA22_genome_delly.tsv")) %>% 
  filter(str_detect(SAMPLE, "1_BM_CD138pos"), 
         SVTYPE=="TRA", 
         (CHROM=="chr14" & CHR2 %in% chr_partners) | (CHROM %in% chr_partners & CHR2=="chr14")) %>%
  mutate(ID = str_remove(SAMPLE, "_1_BM_CD138pos"), 
         CHROM = as.numeric(str_remove(CHROM, "chr")), 
         CHR2 = as.numeric(str_remove(CHR2, "chr"))) %>%
  select(ID, CHROM, POS, REF, CHR2, ENDPOSSV)


# Manta
manta_df <- fread(paste0(filePath, "structural_event/MMRF_CoMMpass_IA22_genome_manta.tsv")) %>% 
  filter(str_detect(SAMPLE, "1_BM_CD138pos"), 
         SVTYPE=="BND", 
         (CHROM=="chr14" & CHR2 %in% chr_partners) | (CHROM %in% chr_partners & CHR2=="chr14")) %>%
  mutate(ID = str_remove(SAMPLE, "_1_BM_CD138pos"), 
         CHROM = as.numeric(str_remove(CHROM, "chr")), 
         CHR2 = as.numeric(str_remove(CHR2, "chr"))) %>%
  select(ID, CHROM, POS, REF, CHR2, ENDPOSSV)

# merge
t_IgH_df <- full_join(delly_df, manta_df, by = c("ID", "CHROM", "POS", "REF", "CHR2", "ENDPOSSV"), relationship = "many-to-many") %>%
  mutate(t_IgH = ifelse(CHROM > CHR2, paste0("t(",CHR2,";",CHROM,")"), paste0("t(",CHROM,";",CHR2,")")), .after = ID)

write_tsv(t_IgH_df, "canonical_t_IgH_NGS.tsv")


