# !/usr/bin/r


# file: 2_mmrf_check_datasets
# aim: checking data consistency
# last update: 06-08-2026


library(data.table)
library(openxlsx)
library(tidyverse)
library(waldo)


# ---- Main ----
# environment setting
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")
outDir <- "/mmrf/sanity_check/"

dir.create(paste0(wd, outDir), recursive = TRUE)


# ---- Loading Datasets ----
# loading original mmrf dataset
mmrf_transcript_based_df <- fread(paste0(wd, "/transcript_based/mmrf_tp53_per_pt.txt")) %>%
  select(PUBLIC_ID, p53_FL, p53_FL_exp, p53_beta, p53_beta_exp, p53_gamma, p53_gamma_exp, delta40_p53_alpha, delta40_p53_alpha_exp, delta133_p53_alpha, 
         delta133_p53_alpha_exp, delta133_p53_beta, delta133_p53_beta_exp, delta133_p53_gamma, delta133_p53_gamma_exp) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(PUBLIC_ID)

mmrf_genomic_based_df <- fread(paste0(wd, "/genomic_based/mmrf_genomic_per_pt_class_10.txt")) %>%
  select(-c(2:42)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(PUBLIC_ID)


# ---- Sanity check: clinical data ----
mmrf_cln_per_pt <- fread("clinical_data/mmrf_cln_per_pt.txt") %>% 
  select(-c(deathdy, lstalive, lvisitdy, lastdy)) %>%
  #mutate(light_chain_type = tolower(light_chain_type)) %>%
  arrange(PUBLIC_ID)

# comparison with mmrf transcript based
compare(mmrf_cln_per_pt, mmrf_transcript_based_df %>% select(all_of(colnames(mmrf_cln_per_pt))))

# comparison with mmrf genomic based
compare(mmrf_cln_per_pt, mmrf_genomic_based_df %>% select(all_of(colnames(mmrf_cln_per_pt))))


# ---- Sanity check: categorization at transcript level ----
# computing median per isoform
median_per_tp53 <- mmrf_trascript_based_categorized %>% 
  summarise(across(-PUBLIC_ID, ~median(.[. > 0])))

# categorization
mmrf_transcript_based_categorized <- mmrf_transcript_based_df %>%
  select(PUBLIC_ID, contains("p53"), -ends_with("exp")) %>%
  mutate(across(.cols = -PUBLIC_ID, 
                .fns = ~ ifelse(.x==0, "absent", ifelse(.x > median_per_tp53[[cur_column()]], "high", "low")),
                .names = "{.col}_exp")) %>%
  arrange(PUBLIC_ID)

# comparison
compare(mmrf_transcript_based_df %>% select(PUBLIC_ID, contains("p53")), 
        mmrf_transcript_based_categorized %>% select(all_of(colnames(mmrf_transcript_based_df %>% select(PUBLIC_ID, contains("p53"))))))


# ---- Sanity check: categorization at genomic level ----
# categorization
mmrf_genomic_based_categorized <- mmrf_genomic_based_df %>%
  select(PUBLIC_ID, c(43:163), -ends_with("alt")) %>%
  mutate(across(.cols = c(2:46), 
                .fns = ~ ifelse(.x >= 2.10, "amp", ifelse(.x <= 1.90, "del", "normal")),
                .names = "{.col}_alt")) %>%
  arrange(PUBLIC_ID)

# comparison at chr level
compare(mmrf_genomic_based_categorized %>% select(PUBLIC_ID, starts_with("[0-9]+"), ends_with("alt")), 
        mmrf_genomic_based_df %>% select(PUBLIC_ID, starts_with("[0-9]+"), ends_with("alt"), -ALT))

# comparison at IGH level
compare(mmrf_genomic_based_categorized %>% select(PUBLIC_ID, ends_with("CALL"), t_IgH, IGH_CHROM, IGH_POS, IGH_CHR2, IGH_POS2), 
        mmrf_genomic_based_df %>% select(PUBLIC_ID, ends_with("CALL"), t_IgH, IGH_CHROM, IGH_POS, IGH_CHR2, IGH_POS2))


# comparison at MUT level
compare(mmrf_genomic_based_categorized %>% select(PUBLIC_ID, TP53_CHROM, TP53_POS, ID, REF, MUTECT2, STRELKA2, VARDICT, OCTOPUS, LANCET, EFFECT, HUGO_ID, ENSG_ID, CDS_ANN, PROT_ANN, ALT_DP, REF_DP, VAF), 
        mmrf_genomic_based_df %>% select(PUBLIC_ID, TP53_CHROM, TP53_POS, ID, REF, MUTECT2, STRELKA2, VARDICT, OCTOPUS, LANCET, EFFECT, HUGO_ID, ENSG_ID, CDS_ANN, PROT_ANN, ALT_DP, REF_DP, VAF))


# ---- Categorization ----
# merge all datasets
mmrf_first_line_per_pt <- left_join(mmrf_transcript_based_df, mmrf_genomic_based_df %>% select(PUBLIC_ID, c(43:163)), by = "PUBLIC_ID")

mmrf_first_line_cat_per_pt <- mmrf_first_line_per_pt %>%
  mutate(across(ends_with("exp"), ~ifelse(.=="absent", 1, 0), .names = "{col}_absent")) %>%
  mutate(across(ends_with("exp"), ~ifelse(.=="low", 1, 0), .names = "{col}_low")) %>%
  mutate(across(ends_with("exp"), ~ifelse(.=="high", 1, 0), .names = "{col}_high")) %>%
  mutate(p53FL_p53beta_balanced = ifelse(p53_FL_exp==p53_beta_exp, 1, 0)) %>%
  mutate(across(ends_with("alt"), ~ifelse(.=="amp", 1, 0), .names= "label_{col}_amp")) %>%
  mutate(across(ends_with("alt"), ~ifelse(.=="del", 1, 0), .names= "label_{col}_del")) %>%
  mutate(TP53_MUT = ifelse(is.na(TP53_POS), 0, 1)) %>%
  mutate(label_3_amp = ifelse(label_3p_alt_amp==1 | label_3q_alt_amp==1, 1, 0),
         label_5_amp = ifelse(label_5p_alt_amp==1 | label_5q_alt_amp==1, 1, 0),
         label_7_amp = ifelse(label_7p_alt_amp==1 | label_7q_alt_amp==1, 1, 0),
         label_9_amp = ifelse(label_9p_alt_amp==1 | label_9q_alt_amp==1, 1, 0),
         label_11_amp = ifelse(label_11p_alt_amp==1 | label_11q_alt_amp==1, 1, 0),
         label_19_amp = ifelse(label_19p_alt_amp==1 | label_19q_alt_amp==1, 1, 0),
         label_21_amp = ifelse(label_21p_alt_amp==1 | label_21q_alt_amp==1, 1, 0),
         tot_amp = label_3_amp + label_5_amp + label_7_amp + label_9_amp + label_11_amp + label_15q_alt_amp + label_19_amp, label_21_amp,
         hd = ifelse(tot_amp > 2, 1, 0), 
         hd_with_19q = ifelse(hd==1 & label_19q_alt_amp==1, 1, 0), 
         hd_with_21q = ifelse(hd==1 & label_21q_alt_amp==1, 1, 0)) %>%
  arrange(PUBLIC_ID)

write.xlsx(mmrf_first_line_cat_per_pt, "sanity_check/mmrf_first_line_per_pt.xlsx") #will follow manual revision


# ---- Sanity check with db of 160925 ----
db_160925_rev <- read.xlsx("sanity_check/revision/db_clinico_genomico_sopravvivenza_1609257_revisionato_070826.xlsx") %>%
  arrange(PUBLIC_ID) %>%
  
  mutate(across(where(is.numeric), round(., 3)))

mmrf_first_line_cat_per_pt_rev <- read.xlsx("sanity_check/mmrf_first_line_per_pt.xlsx") %>%
  select(all_of(colnames(db_160925_rev))) %>%
  arrange(PUBLIC_ID)

compare(db_160925_rev %>% select(PUBLIC_ID, `1p`), mmrf_first_line_cat_per_pt_rev %>% select(PUBLIC_ID, `1p`), max_diffs = Inf)


# ---- Checking isoforms expressions ----
mmrf_tp53_per_pt <- fread("transcript_based/mmrf_tp53_per_pt.txt") %>%
  select(PUBLIC_ID, ends_with("exp")) %>%
  arrange(PUBLIC_ID)

mmrf_first_line_tp53_per_pt <- mmrf_first_line_cat_per_pt_rev %>%
  select(PUBLIC_ID, ends_with("exp")) %>%
  arrange(PUBLIC_ID)

waldo::compare(mmrf_tp53_per_pt, mmrf_first_line_tp53_per_pt)


# ---- Last sanity check -----
mmrf_complete_df <- read.xlsx("sanity_check/mmrf_first_line_per_pt_complete.xlsx")


# check transcript expression
mmrf_complete_transcript_based <- mmrf_complete_df %>%
  select(PUBLIC_ID, p53_FL, p53_FL_exp, p53_beta, p53_beta_exp, p53_gamma, p53_gamma_exp, delta40_p53_alpha, delta40_p53_alpha_exp, delta133_p53_alpha, 
         delta133_p53_alpha_exp, delta133_p53_beta, delta133_p53_beta_exp, delta133_p53_gamma, delta133_p53_gamma_exp) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(PUBLIC_ID)

compare(mmrf_complete_transcript_based, mmrf_transcript_based_df)

# check if categorization is correct
median_per_tp53 <- mmrf_complete_transcript_based %>% 
  summarise(across(where(is.numeric), ~median(.[. > 0])))


# check genomic expression and categorization
mmrf_complete_genomic_df <- mmrf_complete_df %>%
  select(PUBLIC_ID, c(63:183)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(PUBLIC_ID)

compare(mmrf_complete_genomic_df, mmrf_genomic_based_df, max_diffs = Inf)


# check clinical data
mmrf_complete_clinical_df <- mmrf_complete_df %>%
  select(c(1:48), -starts_with("response"), -resp_group) %>%
  arrange(PUBLIC_ID)

compare(mmrf_complete_clinical_df, mmrf_cln_per_pt, max_diffs = Inf)


# check transcripts multiple categorizations
mmrf_transcript_based_cats <- mmrf_complete_df %>%
  select(PUBLIC_ID, p53_FL, p53_FL_exp, p53_beta, p53_beta_exp, p53_gamma, p53_gamma_exp, delta40_p53_alpha, delta40_p53_alpha_exp, delta133_p53_alpha, 
         delta133_p53_alpha_exp, delta133_p53_beta, delta133_p53_beta_exp, delta133_p53_gamma, delta133_p53_gamma_exp, 
         c(184:217)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(PUBLIC_ID)


# check genomic multiple categorizations
mmrf_genomic_based_cats <- mmrf_complete_df %>%
  select(PUBLIC_ID, c(63:183), c(218:319)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(PUBLIC_ID)
