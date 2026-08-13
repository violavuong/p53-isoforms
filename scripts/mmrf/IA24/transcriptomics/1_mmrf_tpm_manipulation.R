# !/usr/bin/r


# file: 1_mmrf_tpm_manipulation
# aim: clean and manipulate TPM-normalized transcriptomic data of MMRF IA24 (IG filtered) in different lines, specimen
# last update: 12-08-26


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/transcriptomic/")
mmrfDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/transcript_expression/"


# ---- MMRF TPM transcriptomic data ----
mmrf_tpm_data <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA24_salmon_transcriptUnstrandedIgFiltered_tpm.tsv"))
dim(mmrf_tpm_data) #224779 transcripts for 933 pts - same pts in different lines

# check if there are unique ENST
mmrf_tpm_data$Transcript %>% unique() %>% length()

# transpose the dataset so that there is a pts per line per cell type for each row, and all transcripts on the columns
# 806 unique pts, 6 possible therapy lines, 2 possible specimen (BM/PB), all CD138pos
mmrf_tpm_per_pt <- mmrf_tpm_data %>%
  column_to_rownames("Transcript") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  mutate(public_id = tolower(str_extract(ID, "MMRF_[0-9]+")), 
         line = str_extract(ID, "_[0-9]_") %>% str_remove_all("_"), 
         specimen = str_extract(ID, "[BMPB]+_CD138pos") %>% str_remove("_CD138pos"), 
         .after = ID) %>%
  select(-ID) %>%
  arrange(public_id, line)

write.table(mmrf_tpm_per_pt, "/complete_db/mmrf_tpm_per_pt_ia24.txt", sep = "\t", dec = ",", na = "", row.names = FALSE, quote = FALSE)
