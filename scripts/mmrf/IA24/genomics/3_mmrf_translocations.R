# !/usr/bin/r


# file: 3_mmrf_translocations 
# aim: extract translocation data from commpass ia24
# last update: 13-08-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/genomics/")
mmrfDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/"


# ---- Load mmrf ia24 traslocations datasets ----
mmrf_igtx <- fread(paste0(mmrfDir, "/Ig_translocations/MMRF_CoMMpass_IA24_genome_tumor_only_mm_igtx_pairoscope.tsv"))
mmrf_delly <- fread(paste0(mmrfDir, "/structural_events/MMRF_CoMMpass_IA24_genome_delly.tsv"))
mmrf_manta <- fread(paste0(mmrfDir, "/structural_events/MMRF_CoMMpass_IA24_genome_manta.tsv"))


# ---- define igh translocations breakpoints (following uma) ----
IgH_df <- IgH_df %>% filter(POS_IgH > 105586437 - IgHLocus & POS_IgH < 106879844 + IgHLocus)


igh_breaks <- data.frame(partner = c("nsd2", "ccnd3", "myc", "mafa", "ccnd1", "ccnd2", "maf", "mafb"), 
                         chr_partner = c(4, 6, 8, 8, 11, 12, 16, 20), 
                         start_hg38 = c(1870559, 41934933, 127735434, 143419182, 69641143, 4269771, 79202622, 40685848), 
                         end_hg38 = c(1982207, 42050357, 127742951, 143430732, 69654474, 4306933, 79600737, 40689236), 
                         mean_pos = c(1926383, 41992645, 127739193, 143424957, 69647809, 4288352, 79401680, 40687542))



# ---- igtx data manipulation ----
# cleaning
mmrf_igtx_per_pt <- mmrf_igtx %>%
  mutate(public_id = tolower(str_extract(SAMPLE, "MMRF_[0-9]+")), 
         line = str_extract(SAMPLE, "(?<=MMRF_[0-9]{4}_)[0-9]+"), 
         specimen = str_extract(SAMPLE, "(BM|PB)_CD138pos") %>% str_remove("_CD138pos")) %>%
  select(public_id, line, specimen, ends_with("CALL"), ends_with("IGSOURCE")) %>%
  arrange(public_id, line)

write.xlsx(mmrf_igtx_per_pt, "mmrf_igtx_per_pt_ia24.xlsx")

# filtering
mmrf_igtx_first_line <- mmrf_igtx_per_pt %>%
  filter(line==1, specimen=="BM") %>%
  mutate(igtx_nsd2 = ifelse(NSD2_CALL==1 & NSD2_IGSOURCE==1, 1, 0),
         igtx_ccnd3 = ifelse(CCND3_CALL==1 & CCND3_IGSOURCE==1, 1, 0), 
         igtx_myc = ifelse(MYC_CALL==1 & MYC_IGSOURCE==1, 1, 0), 
         igtx_mafa = ifelse(MAFA_CALL==1 & MAFA_IGSOURCE==1, 1, 0),
         igtx_ccnd1 = ifelse(CCND1_CALL==1 & CCND1_IGSOURCE==1, 1, 0), 
         igtx_ccnd2 = ifelse(CCND2_CALL==1 & CCND2_IGSOURCE==1, 1, 0), 
         igtx_maf = ifelse(MAF_CALL==1 & MAF_IGSOURCE==1, 1, 0), 
         igtx_mafb = ifelse(MAFB_CALL==1 & MAFB_IGSOURCE==1, 1, 0)) %>%
  select(public_id, starts_with("igtx"))


# ---- delly data manipulation ----
# cleaning
mmrf_delly_per_pt <- mmrf_delly %>%
  mutate(public_id = tolower(str_extract(SAMPLE, "MMRF_[0-9]+")), 
         line = str_extract(SAMPLE, "(?<=MMRF_[0-9]{4}_)[0-9]+"), 
         specimen = str_extract(SAMPLE, "(BM|PB)_CD138pos") %>% str_remove("_CD138pos"), 
         CHROM = str_remove(CHROM, "chr"), 
         CHR2 = str_remove(CHR2, "chr")) %>%
  filter(SVTYPE=="TRA") %>%
  select(public_id, line, specimen, CHROM, POS, CHR2, ENDPOSSV) %>%
  arrange(public_id, line)

write.xlsx(mmrf_delly_per_pt, "mmrf_delly_per_pt_ia24.xlsx")


# filtering
mmrf_delly_per_pt$POS_partner <- ifelse(mmrf_delly_per_pt$CHROM==14, mmrf_delly_per_pt$ENDPOSSV, mmrf_delly_per_pt$POS)
mmrf_delly_per_pt$CHROM_partner <- ifelse(mmrf_delly_per_pt$CHROM==14, mmrf_delly_per_pt$CHR2, mmrf_delly_per_pt$CHROM)

mmrf_delly_per_pt$POS_igh <- ifelse(mmrf_delly_per_pt$CHROM==14, mmrf_delly_per_pt$POS, mmrf_delly_per_pt$ENDPOSSV)
mmrf_delly_per_pt$CHROM_igh <- ifelse(mmrf_delly_per_pt$CHROM==14, mmrf_delly_per_pt$CHROM, mmrf_delly_per_pt$CHR2)

mmrf_delly_first_line <- mmrf_delly_per_pt %>% 
  filter(line==1, specimen=="BM") %>%
  filter(POS_igh > 106052774 - 1e6 & POS_igh < 107288051 + 1e6, 
         CHROM_partner %in% c(4, 6, 8, 11, 12, 16, 20),
         CHROM_igh==14) %>%
  mutate(delly_nsd2 = ifelse(CHROM_partner==4 & (abs(POS_partner - 1926383) < 5e6), 1, 0), 
         delly_ccnd3 = ifelse(CHROM_partner==6 & (abs(POS_partner - 41992645) < 5e6), 1, 0), 
         delly_myc = ifelse(CHROM_partner==8 & (abs(POS_partner - 127739193) < 5e6), 1, 0), 
         delly_mafa = ifelse(CHROM_partner==8 & (abs(POS_partner - 143424957) < 5e6), 1, 0), 
         delly_ccnd1 = ifelse(CHROM_partner==11 & (abs(POS_partner - 69647809) < 5e6), 1, 0), 
         delly_ccnd2 = ifelse(CHROM_partner==12 & (abs(POS_partner - 4288352) < 5e6), 1, 0), 
         delly_maf = ifelse(CHROM_partner==16 & (abs(POS_partner - 79401680) < 5e6), 1, 0), 
         delly_mafb = ifelse(CHROM_partner==20 & (abs(POS_partner - 40687542) < 5e6), 1, 0)) %>%
  select(public_id, starts_with("delly")) %>%
  group_by(public_id) %>%
  summarise(across(everything(), max))


# ---- Cleaning manta data ----
mmrf_manta_per_pt <- mmrf_manta %>%
  mutate(public_id = tolower(str_extract(SAMPLE, "MMRF_[0-9]+")), 
         line = str_extract(SAMPLE, "(?<=MMRF_[0-9]{4}_)[0-9]+"), 
         specimen = str_extract(SAMPLE, "(BM|PB)_CD138pos") %>% str_remove("_CD138pos"), 
         CHROM = str_remove(CHROM, "chr"), 
         CHR2 = str_remove(CHR2, "chr")) %>%
  filter(SVTYPE=="BND") %>%
  select(public_id, line, specimen, REF, ALT, CHROM, POS, CHR2, ENDPOSSV) %>%
  arrange(public_id, line)

write.xlsx(mmrf_manta_per_pt, "mmrf_manta_per_pt_ia24.xlsx")


# filtering
mmrf_manta_per_pt$POS_partner <- ifelse(mmrf_manta_per_pt$CHROM==14, mmrf_manta_per_pt$ENDPOSSV, mmrf_manta_per_pt$POS)
mmrf_manta_per_pt$CHROM_partner <- ifelse(mmrf_manta_per_pt$CHROM==14, mmrf_manta_per_pt$CHR2, mmrf_manta_per_pt$CHROM)

mmrf_manta_per_pt$POS_igh <- ifelse(mmrf_manta_per_pt$CHROM==14, mmrf_manta_per_pt$POS, mmrf_manta_per_pt$ENDPOSSV)
mmrf_manta_per_pt$CHROM_igh <- ifelse(mmrf_manta_per_pt$CHROM==14, mmrf_manta_per_pt$CHROM, mmrf_manta_per_pt$CHR2)

mmrf_manta_first_line <- mmrf_manta_per_pt %>% 
  filter(line==1, specimen=="BM") %>%
  filter(POS_igh > 106052774 - 1e6 & POS_igh < 107288051 + 1e6, 
         CHROM_partner %in% c(4, 6, 8, 11, 12, 16, 20),
         CHROM_igh==14) %>%
  mutate(manta_nsd2 = ifelse(CHROM_partner==4 & (abs(POS_partner - 1926383) < 5e6), 1, 0), 
         manta_ccnd3 = ifelse(CHROM_partner==6 & (abs(POS_partner - 41992645) < 5e6), 1, 0), 
         manta_myc = ifelse(CHROM_partner==8 & (abs(POS_partner - 127739193) < 5e6), 1, 0), 
         manta_mafa = ifelse(CHROM_partner==8 & (abs(POS_partner - 143424957) < 5e6), 1, 0), 
         manta_ccnd1 = ifelse(CHROM_partner==11 & (abs(POS_partner - 69647809) < 5e6), 1, 0), 
         manta_ccnd2 = ifelse(CHROM_partner==12 & (abs(POS_partner - 4288352) < 5e6), 1, 0), 
         manta_maf = ifelse(CHROM_partner==16 & (abs(POS_partner - 79401680) < 5e6), 1, 0), 
         manta_mafb = ifelse(CHROM_partner==20 & (abs(POS_partner - 40687542) < 5e6), 1, 0)) %>%
  select(public_id, starts_with("manta")) %>%
  group_by(public_id) %>%
  summarise(across(everything(), max))

# ---- Integration ----
mmrf_translocations <- full_join(mmrf_igtx_first_line, mmrf_delly_first_line, by = "public_id") %>%
  full_join(mmrf_manta_first_line, by = "public_id") %>%
  select(public_id, ends_with("nsd2"), ends_with("ccnd3"), ends_with("myc"), ends_with("mafa"), ends_with("ccnd1"), ends_with("ccnd2"), ends_with("maf"), ends_with("mafb")) %>%
  arrange(public_id)

mmrf_translocations[is.na(mmrf_translocations)] <- 0

write.xlsx(mmrf_translocations, "mmrf_translocation_per_pt.xlsx")
