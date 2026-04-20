# !/usr/bin/r


# file: relapse
# aim: extract genomic data for relapse pts
# last update: 16-04-2026


library(data.table)
library(GenomicRanges)
library(plyranges)
library(tidyverse)


source("C:/Users/Dell/Desktop/git_projects/unique-molecular-assay/UMA_lib/fun/Popeye2.R")
source("C:/Users/Dell/Desktop/git_projects/TP53/scripts/fun/utils.R")


# ---- Main ----
# environment setting
wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - mmrf/")
libDir <- "C:/Users/Dell/Desktop/git_projects/unique-molecular-assay/UMA_lib/"
mmrfDir <- "C:/Users/Dell/Desktop/"


# cohort 
mmrf_tp53_relapsed_per_pt <- fread("transcript_based/relapse/mmrf_tp53_relapsed_per_pt_line.txt")


# ---- GATK Genome Broad Copy Number ----
# computing abs CN value
seg_CN <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_genome_gatk_cna.seg")) %>%
  mutate(SPECTRUM_SEQ = str_remove(SAMPLE, "_BM_CD138pos"), 
         CN = 2^(Segment_Mean + 1), 
         Chromosome = str_remove(Chromosome, "chr")) %>%
  #filter(SPECTRUM_SEQ %in% mmrf_tp53_relapsed_per_pt$SPECTRUM_SEQ) %>%
  rename_with(~"ID", SAMPLE)

# assigning chrarm with Popeye/manipulating
chrarm_seg_CN_df <- Popeye2(seg_CN, "hg38", removeXY = FALSE)
chrarm_seg_CN_df$ID <- str_remove(chrarm_seg_CN_df$ID, "_[1-9]_BM_CD138pos")
names(chrarm_seg_CN_df)[names(chrarm_seg_CN_df)=="ID"] <- "PUBLIC_ID"

# computing the weighted mean CN value
weighted_CN_per_pt <- chrarm_seg_CN_df %>%
  group_by(PUBLIC_ID, SPECTRUM_SEQ, chr, chrarm) %>%
  summarise(weighted_mean_CN = sum(CN * width) / sum(width), 
            start = min(start), 
            end = max(end), 
            width = end - start,
            probes = sum(Num_Probes), 
            .drop = "group") %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, chr, chrarm, start, end, width, probes, weighted_mean_CN)
weighted_CN_per_pt$chrarm <- factor(weighted_CN_per_pt$chrarm, levels = unique(weighted_CN_per_pt$chrarm))


# ---- Non-synonymous TP53 Mutations ----
NS_TP53_mut <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_combined_vcfmerger2_All_Canonical_NS_Variants.tsv")) %>%
  mutate(SPECTRUM_SEQ = str_remove(SAMPLE, "_BM_CD138pos"), 
         PUBLIC_ID = str_remove(SAMPLE, "_[1-9]_BM_CD138pos"), 
         VAF = `GEN[1].AD[1]`/(`GEN[1].AD[1]` + `GEN[1].AD[0]`)) %>%
  filter(SPECTRUM_SEQ %in% mmrf_tp53_relapsed_per_pt$SPECTRUM_SEQ, `ANN[*].GENE`=="TP53") %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, TP53_POS = POS, ID, REF, ALT, MUTECT2, STRELKA2, VARDICT, OCTOPUS, LANCET, 
         EFFECT = `ANN[*].EFFECT`, HUGO_ID = `ANN[*].GENE`, ENSG_ID = `ANN[*].GENEID`, CDS_ANN = `ANN[*].HGVS_C`, PROT_ANN = `ANN[*].HGVS_P`, 
         ALT_DP = `GEN[1].AD[1]`, REF_DP = `GEN[1].AD[0]`, VAF)

write_tsv(NS_tp53_per_pt, "NS_TP53_per_pt.txt")


# ---- FISH translocations ----
FISH_per_pt <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_genome_tumor_only_mm_igtx_pairoscope.tsv")) %>%
  mutate(SPECTRUM_SEQ = str_remove(SAMPLE, "_BM_CD138pos"), 
         PUBLIC_ID = str_remove(SAMPLE, "_[1-9]_BM_CD138pos")) %>%
  filter(SPECTRUM_SEQ %in% mmrf_tp53_relapsed_per_pt$SPECTRUM_SEQ) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, ends_with("CALL"))

write_tsv(fish_per_pt, "canonical_t_IgH_FISH_per_pt.txt")


# ---- NGS translocations ----
partners <- c("4", "6", "8", "11", "16", "20") # added MYC

# DELLY
DELLY_df <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_genome_delly.tsv")) %>% 
  mutate(SPECTRUM_SEQ = str_remove(SAMPLE, "_BM_CD138pos"), 
         PUBLIC_ID = str_remove(SAMPLE, "_[1-9]_BM_CD138pos"), 
         CHROM = str_remove(CHROM, "chr"), 
         CHR2 = str_remove(CHR2, "chr")) %>%
  filter(SPECTRUM_SEQ %in% mmrf_tp53_relapsed_per_pt$SPECTRUM_SEQ, SVTYPE=="TRA", (CHROM=="14" & CHR2 %in% partners) | (CHROM %in% partners & CHR2=="14")) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, CHROM, POS, CHR2, POS2 = ENDPOSSV)

# Manta
Manta_df <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_genome_manta.tsv")) %>% 
  mutate(SPECTRUM_SEQ = str_remove(SAMPLE, "_BM_CD138pos"), 
         PUBLIC_ID = str_remove(SAMPLE, "_[1-9]_BM_CD138pos"), 
         CHROM = str_remove(CHROM, "chr"), 
         CHR2 = str_remove(CHR2, "chr")) %>%
  filter(SPECTRUM_SEQ %in% mmrf_tp53_relapsed_per_pt$SPECTRUM_SEQ, SVTYPE=="BND", (CHROM=="14" & CHR2 %in% partners) | (CHROM %in% partners & CHR2=="14")) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, CHROM, POS, CHR2, POS2 = ENDPOSSV)


# NGS IGH 
IGH <- full_join(DELLY_df, Manta_df, by = c("PUBLIC_ID", "SPECTRUM_SEQ", "CHROM", "POS", "CHR2", "POS2"), relationship = "many-to-many") %>%
  mutate(CHROM = as.numeric(CHROM), CHR2 = as.numeric(CHR2), 
         t_IGH = ifelse(CHROM > CHR2, paste0("t(",CHR2,";",CHROM,")"), paste0("t(",CHROM,";",CHR2,")")), .after = SPECTRUM_SEQ)

# swapping chrom and pos
is_chr14 <- IGH[["CHROM"]]!= "14"
IGH[is_chr14, c("CHROM", "CHR2")] <- IGH[is_chr14, c("CHR2", "CHROM")]
IGH[is_chr14, c("POS", "POS2")] <- IGH[is_chr14, c("POS2", "POS")]


write_tsv(t_IgH_per_pt, "canonical_t_IgH_NGS_per_pt.txt")


# ---- Classification ----
# transposing weighted mean CN dataframe 
weighted_CN_per_chrarm_relapsed <- pivot_wider(weighted_CN_per_pt, id_cols = c(PUBLIC_ID, SPECTRUM_SEQ), names_from = chrarm, values_from = weighted_mean_CN) 



weighted_CN_class_10 <- classifyGenomicArms(weighted_CN_per_chrarm, 0.10, chr_arms) 
mmrf_genomic_per_pt_class_10 <- mergeGenomicData(weighted_CN_class_10, mmrf_cln_per_pt, fish_per_pt, t_IgH_per_pt, NS_tp53_per_pt)
write_tsv(mmrf_genomic_per_pt_class_10, "mmrf_genomic_per_pt_class_10.txt")
