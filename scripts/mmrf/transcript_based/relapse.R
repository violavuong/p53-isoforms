# !/usr/bin/r


# file: relapse
# aim: check for those pts that changed line (either progression or relapse). first look for transcriptome, then clinical, then genomics
# last update: 17-04-2026


library(data.table)
library(openxlsx)
library(tidyverse)


# ---- Main ----
# environment setting
wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - mmrf/")
outDir <- "/relapse/"

dir.create(paste0(wd, "/transcript_based/", outDir), recursive = TRUE)

# cohort input
mmrf_tp53_cohort <- fread("transcript_based/mmrf_tp53_per_pt.txt")
length(mmrf_tp53_cohort$PUBLIC_ID) # there are 659 pts in total

IDs <- mmrf_tp53_cohort$PUBLIC_ID

# TP53 isoforms ENSEMBL IDs
mmrf_tp53 <- read.xlsx("mmrf_tp53.xlsx", sheet = "isoforms")


# ---- Transcript-based ----
# starting from normalized count per million data
mmrf_cpm <- fread("transcript_based/mmrf_cpm.txt")

# which ID have more than 1 line?
mmrf_tp53_cpm_lines <- mmrf_cpm %>% 
  select(transcript, starts_with(IDs)) %>%
  filter(transcript %in% mmrf_tp53$`ENSEMBL.mRNA`)

# separating lines
first_line <- grep("_1_BM_CD138pos", colnames(mmrf_tp53_cpm_lines), value = TRUE)
FU_lines <- grep("_[2-9]_BM_CD138pos", colnames(mmrf_tp53_cpm_lines), value = TRUE)

# keeping only relapsed/follow up
mmrf_tp53_cpm_relapsed <- as.data.frame(t(mmrf_tp53_cpm_lines %>%
                                       column_to_rownames(var = "transcript"))) %>%
  rownames_to_column("PUBLIC_ID") %>%
  mutate(line = as.numeric(str_match(PUBLIC_ID, "^MMRF_\\d+_(\\d+)_BM_CD138pos")[, 2]), 
         PUBLIC_ID = str_remove(PUBLIC_ID, "_[\\d]_BM_CD138pos"), 
         SPECTRUM_SEQ = paste0(PUBLIC_ID, "_", line)) %>%
  filter(line != 1) %>%
  arrange(PUBLIC_ID) %>%
  relocate(c(SPECTRUM_SEQ, line), .after = PUBLIC_ID)

relapsed_pts <- unique(mmrf_tp53_cpm_relapsed$PUBLIC_ID)

write_tsv(mmrf_tp53_cpm_relapsed, paste0(wd, "/transcript_based/", outDir, "mmrf_tp53_relapsed.txt"))


# ---- Clinical-based ----
mmrfDir <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - CoMMpass_IA22_FlatFiles/"

# biochemical data of pts collected during visits
md_visit <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_PER_PATIENT_VISIT.tsv")) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, ECOG_PERFORMANCEST, 
         D_LAB_chem_albumin, D_LAB_chem_calcium, D_LAB_chem_creatinine, D_LAB_cbc_hemoglobin, D_LAB_chem_ldh, LDH_level, 
         D_LAB_serum_m_protein, D_LAB_urine_24hr_m_protein, D_LAB_cbc_platelet, D_LAB_serum_beta2_microglobulin, D_LAB_serum_c_reactive_protein, 
         D_IM_lightc, D_IM_LIGHT_CHAIN_BY_FLOW, 
         D_IM_CD138_PC_PERCENT)

# pts treatment and therapy data  
md_trt <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_STAND_ALONE_TRTRESP.tsv")) %>%
  mutate(SPECTRUM_SEQ = paste0(PUBLIC_ID, "_", line), .after = PUBLIC_ID) %>%
  mutate(asct = ifelse(str_detect(bmtx_type, "Autologous"), 1, 0), 
         n_asct = ifelse(asct==0, 0, bmtx_seq), 
         PI = ifelse(str_detect(thercat, "Proteasome"), 1, 0), 
         mAb = ifelse(str_detect(thercat, "CD38"), 1, 0), 
         IMID = ifelse(str_detect(thercat, "IMID"), 1, 0), 
         PI_IMID = ifelse(PI==1 & IMID==1, 1, 0)) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, therstdy, therendy, thername, thercat, brespdy, bestresp, bestrespsh, asct, n_asct, PI, mAb, IMID, PI_IMID)

# pts treatment regimen
line_dict <- c("Second" = 2, "Third" = 3, "Fourth" = 4, "Fifth" = 5)

md_trt_reg <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_STAND_ALONE_TREATMENT_REGIMEN.tsv")) %>%
  mutate(line_chr = str_extract(MMTX_TYPE, paste(names(line_dict), collapse = "|")), 
         line = line_dict[line_chr], 
         SPECTRUM_SEQ = paste0(PUBLIC_ID, "_", line), 
         maintenance = ifelse(MMTX_ISTHISLINEOFT=="Maintenance", 1, 0), 
         maint_lena = ifelse(MMTX_ISTHISLINEOFT=="Maintenance" & MMTX_THERAPY=="Lenalidomide", 1, 0), 
         consolidation = ifelse(MMTX_ISTHISLINEOFT=="Consolidation", 1, 0)) %>%
  group_by(PUBLIC_ID, SPECTRUM_SEQ) %>%
  summarise(MMTX_REASONFORTREA = paste(unique(MMTX_REASONFORTREA[!is.na(MMTX_REASONFORTREA) & MMTX_REASONFORTREA != ""]), collapse = "|"), 
            MMTX_ISTHISLINEOFT = paste(unique(MMTX_ISTHISLINEOFT[!is.na(MMTX_ISTHISLINEOFT) & MMTX_ISTHISLINEOFT != ""]), collapse = "|"), 
            maintenance = maintenance[which.max(maintenance)], 
            maint_lena = maint_lena[which.max(maint_lena)], 
            consolidation = consolidation[which.max(consolidation)], 
            .groups = "drop")

# pts survival data  
md_surv <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_STAND_ALONE_SURVIVAL.tsv")) %>%
  select(PUBLIC_ID,
         start_dy_line_2 = linesdy2, end_dy_line_2 = lineedy2, best_respdy_line_2 = bstdy2, 
         start_dy_line_3 = linesdy3, end_dy_line_3 = lineedy3, best_respdy_line_3 = bstdy3, 
         start_dy_line_4 = linesdy4, end_dy_line_4 = lineedy4, best_respdy_line_4 = bstdy4, 
         start_dy_line_5 = linesdy5, end_dy_line_5 = lineedy5, best_respdy_line_5 = bstdy5, 
         PFS_dy_line_2 = pfsdy2, PFS_event_line_2 = censpfs2, PFS_time_line_2 = ttcpfs2,
         PFS_dy_line_3 = pfsdy3, PFS_event_line_3 = censpfs3, PFS_time_line_3 = ttcpfs3,
         OS_date = oscdy, OS_event = censos, OS_time = ttcos, 
         deathdy, lstalive, lvisitdy, lastdy) %>%
  pivot_longer(cols = matches("_line_\\d+$"), names_to = c(".value", "line"), names_pattern = "(.*)_line_(\\d+)$") %>%
  mutate(SPECTRUM_SEQ = paste0(PUBLIC_ID, "_", line)) %>%
  select(PUBLIC_ID, SPECTRUM_SEQ, start_dy, end_dy, best_respdy, PFS_dy, PFS_event, PFS_time, OS_date, OS_event, OS_time, deathdy, lstalive, lvisitdy, lastdy)

# demographic data 
md_pt <- fread(paste0(mmrfDir, "MMRF_CoMMpass_IA22_PER_PATIENT.tsv")) %>%
  select(PUBLIC_ID, D_PT_age, D_PT_gender, D_PT_iss)


# ---- Merging ----
mmrf_tp53_relapsed_per_pt <- left_join(mmrf_tp53_cpm_relapsed, md_pt, by = "PUBLIC_ID") %>%
  left_join(md_visit, by = c("PUBLIC_ID", "SPECTRUM_SEQ")) %>%
  left_join(md_trt, by = c("PUBLIC_ID", "SPECTRUM_SEQ")) %>%
  left_join(md_trt_reg, by = c("PUBLIC_ID", "SPECTRUM_SEQ")) %>%
  left_join(md_surv, by = c("PUBLIC_ID", "SPECTRUM_SEQ")) %>%
  relocate(starts_with("ENST"), .after = lastdy)

write_tsv(mmrf_tp53_relapsed_per_pt, paste0(wd, "/transcript_based/", outDir, "mmrf_tp53_relapsed_per_pt_line.txt"))


# ---- Summarizing relapsed pts ----
# summarizing relapsed pts
relapsed_pts <- mmrf_tp53_relapsed_per_pt %>%
  distinct(PUBLIC_ID, line) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = line, values_from = present)

