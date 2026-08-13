# !/usr/bin/r


# file: 4_mmrf_mutations
# aim: extract and manipulated mmrf ia24 mutations data
# last update: 13-08-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/genomics/")
mmrfDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/mutations/"


# ---- Loading datasets ----
mmrf_ns_mutations <- fread(paste0(mmrfDir, "/MMRF_CoMMpass_IA24_combined_vcfmerger2_All_Canonical_NS_Variants.tsv"))

# ---- data manipulation ----
mmrf_ns <- mmrf_ns_mutations %>%
  mutate(public_id = tolower(str_extract(SAMPLE, "MMRF_[0-9]+")), 
         line = str_extract(SAMPLE, "(?<=MMRF_[0-9]{4}_)[0-9]+"), 
         specimen = str_extract(SAMPLE, "(BM|PB)_CD138pos") %>% str_remove("_CD138pos"), 
         CHROM = str_remove(CHROM, "chr"), 
         VAF = `GEN[1].AD[1]`/(`GEN[1].AD[1]` + `GEN[1].AD[0]`)) %>%
  #filter(`ANN[*].BIOTYPE`=="protein_coding") %>%
  select(public_id, line, specimen, CHROM, POS, REF, ALT, COSMIC, COSMIC_NC, DB, GNOMAD_EXOME, GNOMAD_GENOME, CLINVAR, MUTECT2, STRELKA2, VARDICT, OCTOPUS, LANCET, 
         `ANN[*].EFFECT`, `ANN[*].IMPACT`, `ANN[*].GENE`, `ANN[*].BIOTYPE`, `ANN[*].HGVS_C`, `ANN[*].HGVS_P`, VAF) %>%
  arrange(public_id, line)

mmrf_ns$`ANN[*].GENE` <- factor(mmrf_ns$`ANN[*].GENE`, levels = unique(sort(mmrf_ns$`ANN[*].GENE`)))

write.xlsx(mmrf_ns, "complete_db/mmrf_ns_mut_per_pt_ia24.xlsx")

# create counter dataset
mmrf_ns_mutations_counter_per_pt <- mmrf_ns %>%
  mutate(presence = 1) %>%
  select(all_of(c("public_id", "line", "specimen")), `ANN[*].GENE`, presence) %>%
  distinct() %>%
  pivot_wider(names_from = `ANN[*].GENE`, values_from = presence, values_fill = 0, names_sort = TRUE)

write.table(mmrf_ns_mutations_counter_per_pt, "complete_db/mmrf_ns_mut_counter_ia24.txt", sep = "\t", dec = ",", na = "", row.names = FALSE, quote = FALSE)


# ---- filtering ----
mmrf_ns_first_line <- mmrf_ns_mutations_counter_per_pt %>%
  filter(line==1, specimen=="BM") %>%
  select(-line, -specimen)

write.table(mmrf_ns_first_line, "mmrf_ns_mut_per_pt.txt", sep = "\t", dec = ",", na = "", row.names = FALSE, quote = FALSE)
