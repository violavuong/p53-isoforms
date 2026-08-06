# !/usr/bin/r

# file: 1_mmrf_extraction
# description: data extraction and handling of compelling CoMMpass genomic data
# last update: 06-08-2026

library(data.table)
library(GenomicRanges)
library(plyranges)
library(rlist)
library(tidyverse)

source("C:/Users/violameixian.vuong2/Desktop/git_projects/UMA/UMA_lib/fun/Popeye2.R")
source("C:/Users/violameixian.vuong2/Desktop/git_projects/TP53/scripts/fun/utils.R")

# ---- Main ----
# working directories
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/genomic_based/")
filePath <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA22/"
libDir <- "C:/Users/violameixian.vuong2/Desktop/git_projects/UMA/UMA_lib/"


# ---- GATK genome copy number ----
gene_CN_df <- fread(paste0(filePath, "copy_number/MMRF_CoMMpass_IA22_genome_gatk_cna.seg")) %>%
  filter(str_detect(SAMPLE, "1_BM_CD138pos")) %>%
  mutate(CN = 2^(Segment_Mean + 1), 
         Chromosome = str_remove(Chromosome, "chr")) %>%
  rename_with(~"ID", SAMPLE)

chrarm_CN_df <- Popeye2(gene_CN_df, "hg38", removeXY = FALSE)
chrarm_CN_df$ID <- str_remove(chrarm_CN_df$ID, "_1_BM_CD138pos")
names(chrarm_CN_df)[names(chrarm_CN_df)=="ID"] <- "PUBLIC_ID"

# computing the weighted mean CN value: 877 unique pts
weighted_CN_per_pt <- chrarm_CN_df %>%
  group_by(PUBLIC_ID, chr, chrarm) %>%
  summarise(weighted_mean_CN = sum(CN * width) / sum(width), 
            start = min(start), 
            end = max(end), 
            width = end - start,
            probes = sum(Num_Probes)) %>%
  select(PUBLIC_ID, chr, chrarm, start, end, width, probes, weighted_mean_CN)
weighted_CN_per_pt$chrarm <- factor(weighted_CN_per_pt$chrarm, levels = unique(weighted_CN_per_pt$chrarm))

write_tsv(weighted_CN_per_pt, "genomic_based/genomic_alt/WGS_weighted_CN_per_pt.txt")


# ---- Non-synonymous SNVs ----
# 57 obs for 48 unique pts
tp53_nonsyn_snv_tmp_df <- fread(paste0(filePath, "somatic_mutations/MMRF_CoMMpass_IA22_combined_vcfmerger2_All_Canonical_NS_Variants.tsv")) %>%
  filter(str_detect(SAMPLE, "1_BM_CD138pos"), `ANN[*].GENE`=="TP53") %>%
  mutate(PUBLIC_ID = str_remove(SAMPLE, "_1_BM_CD138pos"), 
         VAF = `GEN[1].AD[1]`/(`GEN[1].AD[1]` + `GEN[1].AD[0]`)) %>%
  select(PUBLIC_ID, TP53_CHROM = CHROM, TP53_POS = POS, ID, REF, ALT, MUTECT2, STRELKA2, VARDICT, OCTOPUS, LANCET, 
         EFFECT = `ANN[*].EFFECT`, HUGO_ID = `ANN[*].GENE`, ENSG_ID = `ANN[*].GENEID`, CDS_ANN = `ANN[*].HGVS_C`, PROT_ANN = `ANN[*].HGVS_P`, 
         ALT_DP = `GEN[1].AD[1]`, REF_DP = `GEN[1].AD[0]`, VAF)

tp53_dups <- tp53_nonsyn_snv_tmp_df[duplicated(tp53_nonsyn_snv_tmp_df$PUBLIC_ID) | duplicated(tp53_nonsyn_snv_tmp_df$PUBLIC_ID, fromLast = TRUE), ] %>%
  group_by(PUBLIC_ID, TP53_CHROM, HUGO_ID, ENSG_ID) %>%
  summarise(across(everything(), ~ paste(., collapse = "/")), .groups = 'drop')

NS_tp53_per_pt <- rbind(tp53_nonsyn_snv_tmp_df %>% filter(!PUBLIC_ID %in% tp53_dups$PUBLIC_ID), tp53_dups)
NS_tp53_per_pt$TP53_CHROM <- as.numeric(str_remove(NS_tp53_per_pt$TP53_CHROM, "chr"))
 
write_tsv(NS_tp53_per_pt, "genomic_based/genomic_alt/NS_TP53_per_pt.txt")


# ---- Translocations ----
chr_partners <- c("chr4", "chr6", "chr11", "chr16", "chr20")

# FISH: 908 unique obs for 908 unique pts
fish_per_pt <- fread(paste0(filePath, "seqFISH/MMRF_CoMMpass_IA22_genome_tumor_only_mm_igtx_pairoscope.tsv")) %>%
  filter(str_detect(SAMPLE, "1_BM_CD138pos")) %>%
  mutate(PUBLIC_ID = str_remove(SAMPLE, "_1_BM_CD138pos")) %>%
  select(PUBLIC_ID, ends_with("CALL"))

write_tsv(fish_per_pt, "genomic_based/genomic_alt/canonical_t_IgH_FISH_per_pt.txt")

# NGS: 1461 obs for 315 unique pts
violameixian.vuong2y_df <- fread(paste0(filePath, "structural_event/MMRF_CoMMpass_IA22_genome_delly.tsv")) %>% 
  filter(str_detect(SAMPLE, "1_BM_CD138pos"), SVTYPE=="TRA", (CHROM=="chr14" & CHR2 %in% chr_partners) | (CHROM %in% chr_partners & CHR2=="chr14")) %>%
  mutate(PUBLIC_ID = str_remove(SAMPLE, "_1_BM_CD138pos"), 
         CHROM = as.numeric(str_remove(CHROM, "chr")), 
         CHR2 = as.numeric(str_remove(CHR2, "chr"))) %>%
  select(PUBLIC_ID, CHROM, POS, CHR2, POS2 = ENDPOSSV)

manta_df <- fread(paste0(filePath, "structural_event/MMRF_CoMMpass_IA22_genome_manta.tsv")) %>% 
  filter(str_detect(SAMPLE, "1_BM_CD138pos"), SVTYPE=="BND", (CHROM=="chr14" & CHR2 %in% chr_partners) | (CHROM %in% chr_partners & CHR2=="chr14")) %>%
  mutate(PUBLIC_ID = str_remove(SAMPLE, "_1_BM_CD138pos"), 
         CHROM = as.numeric(str_remove(CHROM, "chr")), 
         CHR2 = as.numeric(str_remove(CHR2, "chr"))) %>%
  select(PUBLIC_ID, CHROM, POS, CHR2, POS2 = ENDPOSSV)

t_IgH_tmp <- full_join(delly_df, manta_df, by = c("PUBLIC_ID", "CHROM", "POS", "CHR2", "POS2"), relationship = "many-to-many") %>%
  mutate(t_IgH = ifelse(CHROM > CHR2, paste0("t(",CHR2,";",CHROM,")"), paste0("t(",CHROM,";",CHR2,")")), .after = PUBLIC_ID)

# arranging to check for duplicates
is_IgH_chr <- t_IgH_tmp[["CHROM"]]!= "14"
t_IgH_tmp[is_IgH_chr, c("CHROM", "CHR2")] <- t_IgH_tmp[is_IgH_chr, c("CHR2", "CHROM")] #swapping chromosomes
t_IgH_tmp[is_IgH_chr, c("POS", "POS2")] <- t_IgH_tmp[is_IgH_chr, c("POS2", "POS")] #swapping positions

# cleaning 
t_IgH_tmp <- t_IgH_tmp %>%
  filter(POS >= 105586437 & POS <= 106879844) %>% #only falling inside IgH
  distinct(PUBLIC_ID, t_IgH, CHROM, POS, CHR2, .keep_all = TRUE) #removing equal dups

t_IgH_dups <- t_IgH_tmp[duplicated(t_IgH_tmp$PUBLIC_ID) | duplicated(t_IgH_tmp$PUBLIC_ID, fromLast = TRUE), ] %>%
  group_by(PUBLIC_ID) %>%
  summarise(across(everything(), ~ paste(., collapse = "/")), .groups = 'drop') 

# binding
t_IgH_per_pt <- rbind(t_IgH_tmp %>% filter(!PUBLIC_ID %in% t_IgH_dups$PUBLIC_ID), t_IgH_dups) %>%
  rename_with(~c("IGH_CHROM", "IGH_POS", "IGH_CHR2", "IGH_POS2"), c(3,4,5,6))

write_tsv(t_IgH_per_pt, "genomic_based/genomic_alt/canonical_t_IgH_NGS_per_pt.txt")


# ---- Classification ----
# transposing weighted mean CN dataframe 
weighted_CN_per_chrarm <- dcast(setDT(weighted_CN_per_pt), PUBLIC_ID ~ chrarm, value.var = "weighted_mean_CN")

# applying different cut-offs
chr_arms <- colnames(weighted_CN_per_chrarm)[2:46]

weighted_CN_class_10 <- classifyGenomicArms(weighted_CN_per_chrarm, 0.10, chr_arms) 
weighted_CN_class_20 <- classifyGenomicArms(weighted_CN_per_chrarm, 0.20, chr_arms) 
weighted_CN_class_50 <- classifyGenomicArms(weighted_CN_per_chrarm, 0.50, chr_arms) 
weighted_CN_class_80 <- classifyGenomicArms(weighted_CN_per_chrarm, 0.80, chr_arms) 


# ---- Harmonization ----
# 754 pts with clinical and transcriptomic data attached
mmrf_cln_per_pt <- fread("../clinical_data/mmrf_cln_per_pt.txt")
mmrf_cln_per_pt <- mmrf_cln_per_pt %>% filter(mmrf_cln_per_pt$PUBLIC_ID %in% weighted_CN_per_pt$PUBLIC_ID) #filtering out pts that do not have CN

# update
write_tsv(mmrf_cln_per_pt, "../clinical_data/mmrf_cln_per_pt.txt")

# merging
mmrf_genomic_per_pt_class_10 <- mergeGenomicData(weighted_CN_class_10, mmrf_cln_per_pt, fish_per_pt, t_IgH_per_pt, NS_tp53_per_pt)
write_tsv(mmrf_genomic_per_pt_class_10, "mmrf_genomic_per_pt_class_10.txt")

mmrf_genomic_per_pt_class_20 <- mergeGenomicData(weighted_CN_class_20, mmrf_cln_per_pt, fish_per_pt, t_IgH_per_pt, NS_tp53_per_pt)
write_tsv(mmrf_genomic_per_pt_class_20, "mmrf_genomic_per_pt_class_20.txt")

mmrf_genomic_per_pt_class_50 <- mergeGenomicData(weighted_CN_class_50, mmrf_cln_per_pt, fish_per_pt, t_IgH_per_pt, NS_tp53_per_pt)
write_tsv(mmrf_genomic_per_pt_class_50, "mmrf_genomic_per_pt_class_50.txt")

mmrf_genomic_per_pt_class_80 <- mergeGenomicData(weighted_CN_class_80, mmrf_cln_per_pt, fish_per_pt, t_IgH_per_pt, NS_tp53_per_pt)
write_tsv(mmrf_genomic_per_pt_class_80, "mmrf_genomic_per_pt_class_80.txt")
