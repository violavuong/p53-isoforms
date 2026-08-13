#!/usr/bin/r

# file: 1_mmrf_normalization
# aim: investigating mmrf transcript-based data, normalizing and merging clinical data
# last update: 11-08-2026

# installing RNA-seq libraries
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install(c("edgeR", "Glimma"))


# ---- Main ----
library(data.table)
library(edgeR)
library(Glimma)
library(limma)
library(SummarizedExperiment)
library(tidyverse)


#pathDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA22/"
pathDir <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA24/"
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")

# mmrf input data - normalized transcript-per-million
#mmrf_tpm <- fread(paste0(pathDir, "expression_estimates_transcript_based/MMRF_CoMMpass_IA22_salmon_transcriptUnstrandedIgFiltered_tpm.tsv"))
mmrf_tpm <- fread(paste0(pathDir, "transcript_expression/MMRF_CoMMpass_IA24_salmon_transcriptUnstrandedIgFiltered_tpm.tsv"))

# mmrf clinical data
## md_pt: pts general data
md_pt <- fread(paste0(pathDir, "clinical_flat_files/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_PER_PATIENT.tsv"))

## md_visit: biochemical data of pts collected during visits
md_visit <- fread(paste0(pathDir, "clinical_flat_files/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_PER_PATIENT_VISIT.tsv"))

## md_trt: pts treatment and therapy data
md_trt <- fread(paste0(pathDir, "clinical_flat_files/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_STAND_ALONE_TRTRESP.tsv"))

## md_trt_reg: pts treatment regimen
md_trt_reg <- fread(paste0(pathDir, "clinical_flat_files/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_STAND_ALONE_TREATMENT_REGIMEN.tsv"))


## md_surv: pts survival data
md_surv <- fread(paste0(pathDir, "clinical_flat_files/CoMMpass_IA22_FlatFiles/MMRF_CoMMpass_IA22_STAND_ALONE_SURVIVAL.tsv"))


# ----- MMRF data exploration ---- 
str(mmrf_tpm) #transcripts per million are decimals, transcripts are chars
colnames(mmrf_tpm) 
dim(mmrf_tpm) #224779 transcripts for 933 pts - with replicates


length(mmrf_pt$PUBLIC_ID) #1143 pts with clinical data
str(mmrf_visit)
str(mmrf_trt)


# ---- MMRF clinical data manipulation ----
## pts general data
mmrf_pt <- md_pt %>% select(PUBLIC_ID, D_PT_age, D_PT_gender, D_PT_iss)

## pts biochemical data
mmrf_tmp_visit <- md_visit %>% 
  filter(str_detect(SPECTRUM_SEQ, "_1") & VJ_INTERVAL=="Baseline") %>% 
  select(PUBLIC_ID, ECOG_PERFORMANCEST, 
         D_LAB_chem_albumin, D_LAB_chem_calcium, D_LAB_chem_creatinine, D_LAB_cbc_hemoglobin, D_LAB_chem_ldh, LDH_level, 
         D_LAB_serum_m_protein, D_LAB_urine_24hr_m_protein, D_LAB_cbc_platelet, D_LAB_serum_beta2_microglobulin, D_LAB_serum_c_reactive_protein, 
         D_IM_lightc, D_IM_LIGHT_CHAIN_BY_FLOW, 
         D_IM_CD138_PC_PERCENT)
mmrf_visit <- mmrf_tmp_visit[!duplicated(mmrf_tmp_visit$PUBLIC_ID), ]

## pts treatment data
mmrf_trt <- md_trt %>%
  filter(line==1) %>%
  mutate(asct = ifelse(str_detect(bmtx_type, "Autologous"), 1, 0), 
          n_asct = ifelse(asct==0, 0, bmtx_seq)) %>%
  group_by(PUBLIC_ID) %>%
  summarize(therendy_max = max(therendy, na.rm = TRUE),
            thername = thername[which.max(therendy)],
            thercat = thercat[which.max(therendy)],
            resp = unique(bestresp), 
            respsh = unique(bestrespsh),
            asct = asct[which.max(asct)],
            n_asct = n_asct[which.max(n_asct)]) %>%
  mutate(PI = ifelse(str_detect(thercat, "Proteasome"), 1, 0), 
         mAb = ifelse(str_detect(thercat, "CD38"), 1, 0), 
         IMID = ifelse(str_detect(thercat, "IMID"), 1, 0), 
         PI_IMID = ifelse(PI==1 & IMID==1, 1, 0)) 

## pts regimen data
mmrf_reg <- md_trt_reg %>%
  filter(MMTX_TYPE=="First line of therapy") %>%
  mutate(maintenance = ifelse(MMTX_ISTHISLINEOFT=="Maintenance", 1, 0), 
         maint_lena = ifelse(MMTX_ISTHISLINEOFT=="Maintenance" & MMTX_THERAPY=="Lenalidomide", 1, 0), 
         consolidation = ifelse(MMTX_ISTHISLINEOFT=="Consolidation", 1, 0)) %>%
  group_by(PUBLIC_ID) %>%
  summarise(maintenance = maintenance[which.max(maintenance)], 
            maint_lena = maint_lena[which.max(maint_lena)], 
            consolidation = consolidation[which.max(consolidation)])

## pts survival data
mmrf_surv <- md_surv %>%
  select(PUBLIC_ID, start_line_1 = linesdy1, end_line_1 = lineedy1, best_resp_dy_line_1 = bstdy1, 
         PFS_date = pfsdy1, PFS_event = censpfs1, PFS_time = ttcpfs1, 
         OS_date = oscdy, OS_event = censos, OS_time = ttcos, 
         deathdy, lstalive, lvisitdy, lastdy)
  
## merge - 1143 pts have clinical data attached
mmrf_cln <- left_join(mmrf_pt, mmrf_visit, by = "PUBLIC_ID") %>%
  left_join(mmrf_trt, by = "PUBLIC_ID") %>%
  left_join(mmrf_reg, by = "PUBLIC_ID") %>%
  left_join(mmrf_surv, by = "PUBLIC_ID")

write_tsv(mmrf_cln, paste0(wd, "clinical_data/mmrf_cln.txt"))


# ---- MMRF data normalization ----
mmrf_gse <- SummarizedExperiment(assays = list(counts = mmrf_tpm[,-1]))
mmrf_tmm <- calcNormFactors(mmrf_gse, method = "TMM")

# counts-per-million (normalization by sample)
mmrf_cpm <- cpm(mmrf_tmm)
colnames(mmrf_cpm) <- colnames(mmrf_tpm[, -1])
mmrf_cpm <- as.data.frame(mmrf_cpm) %>% mutate(transcript = mmrf_tpm$Transcript)
write_tsv(mmrf_cpm, paste0(wd, "transcript_based/counts_per_million/mmrf_cpm.txt"))

# counts-per-million in logarithmic scale (normalization by sample)
mmrf_lcpm <- cpm(mmrf_tmm, log = TRUE)
colnames(mmrf_lcpm) <- colnames(mmrf_tpm[, -1])
mmrf_lcpm <- as.data.frame(mmrf_lcpm) %>% mutate(transcript = mmrf_tpm$Transcript)
write_tsv(mmrf_lcpm, paste0(wd, "transcript_based/counts_per_million/mmrf_lcpm.txt"))


# ---- MMRF data filtering ----
# removing replicates, keeping only pts with clinical data attached, PUBLIC_ID as a pt-specific ID: 754 pts
mmrf_cpm_per_pt_tmp <- filterMMRF(mmrf_cpm, mmrf_cln) 
mmrf_cpm_per_pt <- mmrf_cpm_per_pt_tmp %>%
  select(order(colnames(mmrf_cpm_per_pt_tmp))) %>%
  mutate(transcript = mmrf_tpm$Transcript, .before = MMRF_1021)
write_tsv(mmrf_cpm_per_pt, paste0(wd, "transcript_based/counts_per_million/mmrf_cpm_per_pt.txt"))

mmrf_lcpm_per_pt_tmp <- filterMMRF(mmrf_lcpm, mmrf_cln) 
mmrf_lcpm_per_pt <- mmrf_lcpm_per_pt_tmp %>%
  select(order(colnames(mmrf_lcpm_per_pt_tmp))) %>%
  mutate(transcript = mmrf_tpm$Transcript, .before = MMRF_1021)
write_tsv(mmrf_lcpm_per_pt, paste0(wd, "transcript_based/counts_per_million/mmrf_lcpm_per_pt.txt"))


# ---- MMRF clinical data filtering ----
# filtering out pts with transcriptomic data but without clinical data
mmrf_cln_per_pt <- mmrf_cln %>% 
  filter(PUBLIC_ID %in% colnames(mmrf_cpm_per_pt)[-1]) %>%
  rename_with(~c("age", "gender", "ISS", "ECOG",
                 "ALB", "Ca", "creatinine", "HE", "LDH", "LDH_level", 
                 "serum_M_prot", "urine_24h_M_prot", "PLT", "serum_B2_mglob", "serum_PCR", 
                 "light_chain", "light_chain_type", 
                 "CD138_perc"), 
              c(2:19)) %>%
  arrange(PUBLIC_ID)

write_tsv(mmrf_cln_per_pt, paste0(wd, "clinical_data/mmrf_cln_per_pt.txt"))


# update: harmonizing every data: pts with clinical, genomic and transcriptomic data: 659 pts
mmrf_cln_per_pt <- fread(paste0(wd, "clinical_data/mmrf_cln_per_pt.txt"))
mmrf_cpm_per_pt <- fread(paste0(wd, "transcript_based/counts_per_million/mmrf_cpm_per_pt.txt"))
mmrf_lcpm_per_pt <- fread(paste0(wd, "transcript_based/counts_per_million/mmrf_lcpm_per_pt.txt"))

mmrf_cpm_per_pt <- mmrf_cpm_per_pt %>% select(transcript, matches(mmrf_cln_per_pt$PUBLIC_ID))
write_tsv(mmrf_cpm_per_pt, paste0(wd, "transcript_based/counts_per_million/mmrf_cpm_per_pt.txt"))

mmrf_lcpm_per_pt <- mmrf_lcpm_per_pt %>% select(transcript, matches(mmrf_cln_per_pt$PUBLIC_ID))
write_tsv(mmrf_lcpm_per_pt, paste0(wd, "transcript_based/counts_per_million/mmrf_lcpm_per_pt.txt"))

