# !/usr/bin/r

# file: 1_mmrf_counts_manipulation
# aim: clean, manipulate mmrf raw counts dataframe for gsea analyses
# last update: 13-08-26


library(data.table)
library(EnsDb.Hsapiens.v86)
library(tidyverse)


# ---- Main ----
# setting env
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/IA24/gene_expression/")
mmrfDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/gene_expression/"


# ---- Cohort ----
mmrf_gene_expression <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA24_salmon_geneUnstranded_counts.tsv")) %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>% 
  mutate(public_id = tolower(str_extract(ID, "MMRF_[0-9]+")), 
       line = str_extract(ID, "(?<=MMRF_[0-9]{4}_)[0-9]+"), 
       specimen = str_extract(ID, "(BM|PB)_CD138pos") %>% str_remove("_CD138pos")) %>%
  relocate(c(public_id, line, specimen), .after = ID) %>%
  select(-ID) %>%
  arrange(public_id, line)

write.table(mmrf_gene_expression, "complete_db/mmrf_gene_expression_per_pt_ia24.txt", sep = "\t", dec = ",", na = "", row.names = FALSE, quote = FALSE)


# ---- Filtering for first line ----
mmrf_gene_exp_first_line <- mmrf_gene_expression %>%
  filter(line==1, specimen=="BM") %>%
  select(-line, -specimen) %>%
  column_to_rownames("public_id") %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("ENSG_ID")


# ---- adding hugo nomenclature ----
# hg38 ensg-hugo dictionary
ensg_hugo_dict_hg38 <- select(EnsDb.Hsapiens.v86, keys = mmrf_gene_exp_first_line$ENSG_ID, keytype = "GENEID", columns = c("GENEID", "SYMBOL")) # 56871 entries

# merge
mmrf_ensg_hugo_gene_exp_first_line <- left_join(mmrf_gene_exp_first_line, ensg_hugo_dict_hg38,  by = c("ENSG_ID"="GENEID")) %>%
  relocate(SYMBOL, .before = ENSG_ID) %>%
  dplyr::rename(HUGO_ID=SYMBOL)

write.table(mmrf_ensg_hugo_gene_exp_first_line, "mmrf_tp53_gene_expression_per_pt.txt", sep = "\t", dec = ",", na = "", row.names = FALSE, quote = FALSE)

