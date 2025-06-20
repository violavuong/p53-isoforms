#!/usr/bin/r

## file: DESeq.R
## last update: 04-06-2025

# ---- Main ----
library(data.table)
library(DESeq2)
library(SummarizedExperiment)
library(tidyverse)


pathDir <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA22/"
setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/MM group - Vuong_Viola_Meixian/TP53/")

# input data
#mmrf_raw_counts <- fread(paste0(pathDir, "expression_estimates_transcript_based/MMRF_CoMMpass_IA22_salmon_transcriptUnstrandedIgFiltered_counts.tsv")) %>%
#  rename_with(~"transcript", "Transcript")
mmrf_cpm_per_pt <- fread("mmrf_cpm_per_pt.txt")
mmrf_cln_per_pt <- fread("mmrf_cln_per_pt.txt")


# ---- data preparation ----
#transcripts <- mmrf_raw_counts$transcript
transcripts <- mmrf_cpm_per_pt$transcript

# deseq requires only integers
#mmrf_raw_counts <- filterMMRF(mmrf_raw_counts, mmrf_cln_per_pt) %>% mutate(transcript = transcripts)
#mmrf_counts <- round(mmrf_raw_counts[,-755]) %>% mutate(transcript = transcripts, .before = MMRF_1021)

mmrf_cpm_counts <- round(mmrf_cpm_per_pt[,-1]) %>% mutate(transcript = transcripts, .before = MMRF_1021)

# checking PUBLIC_ID order 
all(sort(colnames(mmrf_counts)[-1])==mmrf_cln_per_pt$PUBLIC_ID)

# ---- DESeq2 function ----
runDESeq2 <- function(counts, metadata, ft){
  ft <- !!sym(ft)
  
  metadata$ft <- factor(metadata$ft, levels = unique(metadata$ft))
  dds_tmp_ft <- DESeqDataSetFromMatrix(countData = counts[, -1], colData = metadata, design = ~ft)
  
  # removing low-expressed genes
  exp_genes <- rowSums(counts(dds_tmp_ft)) >= 5
  dds_ft <- dds_tmp_ft[exp_genes,]
  
  dds_ft <- DESeq(dds_ft)
  dds_ft_res <- results(dds_ft, contrast = c(paste0(ft), "0", "1"), alpha = 0.05)
  summary(dds_ft_res)
  
  return(dds_ft_res)
}


# ---- DESeq2 analysis: ASCT ----
## ASCT 1/0
#runDESeq2(mmrf_counts, mmrf_cln_per_pt, "asct")
mmrf_cln_per_pt$asct <- factor(mmrf_cln_per_pt$asct, levels = unique(mmrf_cln_per_pt$asct))
dds_tmp_asct <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts[, -1], colData = mmrf_cln_per_pt, design = ~asct)

# removing low-expressed genes
asct_exp_genes <- rowSums(counts(dds_tmp_asct)) >= 5
dds_asct <- dds_tmp_asct[asct_exp_genes,]

dds_asct <- DESeq(dds_asct)
dds_asct_res <- results(dds_asct, contrast = c("asct", "0", "1"), alpha = 0.05)
summary(dds_asct_res)

write.csv(dds_asct_res, "dds_asct.csv")

## number of ASCT (1 vs 2)
# subsetting: taking only pt with asct
mmrf_cln_per_pt_asct <- mmrf_cln_per_pt %>% filter(asct==1)
mmrf_cpm_counts_asct <- mmrf_cpm_counts %>% select(any_of(mmrf_cln_per_pt_asct$PUBLIC_ID))

mmrf_cln_per_pt_asct$n_asct <- factor(ifelse(mmrf_cln_per_pt_asct$n_asct==1, 1, 2), levels = c(1, 2))
dds_tmp_n_asct <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts_asct, colData = mmrf_cln_per_pt_asct, design = ~n_asct)

# removing low-expressed genes
n_asct_exp_genes <- rowSums(counts(dds_tmp_n_asct)) >= 5
dds_n_asct <- dds_tmp_n_asct[n_asct_exp_genes,]

dds_n_asct <- DESeq(dds_n_asct)
dds_n_asct_res <- results(dds_n_asct, contrast = c("n_asct", "1", "2"), alpha = 0.05)
summary(dds_n_asct_res)

write.csv(dds_n_asct_res, "dds_n_asct.csv")

# ---- DESeq2 analysis: maintenance ----
## maintenance 1/0
mmrf_cln_per_pt$maintenance <- factor(mmrf_cln_per_pt$maintenance, levels = unique(mmrf_cln_per_pt$maintenance))
dds_tmp_mant <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts[, -1], colData = mmrf_cln_per_pt, design = ~maintenance)

# removing low-expressed genes
mant_exp_genes <- rowSums(counts(dds_tmp_mant)) >= 5
dds_mant <- dds_tmp_mant[mant_exp_genes,]

dds_mant <- DESeq(dds_mant)
dds_mant_res <- results(dds_mant, contrast = c("maintenance", "0", "1"), alpha = 0.05)
summary(dds_mant_res)

write.csv(dds_mant_res, "dds_mant.csv")

## maintenance lena 1/0
# subsetting: taking only pt with maintenance
mmrf_cln_per_pt_mant <- mmrf_cln_per_pt %>% filter(maintenance==1)
mmrf_cpm_counts_mant <- mmrf_cpm_counts %>% select(any_of(mmrf_cln_per_pt_mant$PUBLIC_ID))

mmrf_cln_per_pt_mant$mant_lena <- factor(mmrf_cln_per_pt_mant$mant_lena, levels = unique(mmrf_cln_per_pt_mant$mant_lena))
dds_tmp_mant_lena <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts_mant, colData = mmrf_cln_per_pt_mant, design = ~mant_lena)

# removing low-expressed genes
mant_lena_exp_genes <- rowSums(counts(dds_tmp_mant_lena)) >= 5
dds_mant_lena <- dds_tmp_mant_lena[mant_lena_exp_genes,]

dds_mant_lena <- DESeq(dds_mant_lena)
dds_mant_lena_res <- results(dds_mant_lena, contrast = c("mant_lena", "0", "1"), alpha = 0.05)
summary(dds_mant_lena_res)

write.csv(dds_mant_lena_res, "dds_mant_lena.csv")

# ---- DESeq2 analysis: consolidation ----
mmrf_cln_per_pt$consolidation <- factor(mmrf_cln_per_pt$consolidation, levels = unique(mmrf_cln_per_pt$consolidation))
dds_tmp_cons <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts[, -1], colData = mmrf_cln_per_pt, design = ~consolidation)

# removing low-expressed genes
cons_exp_genes <- rowSums(counts(dds_tmp_cons)) >= 5
dds_cons <- dds_tmp_cons[cons_exp_genes,]

dds_cons <- DESeq(dds_cons)
dds_cons_res <- results(dds_cons, contrast = c("consolidation", "0", "1"), alpha = 0.05)
summary(dds_cons_res)

write.csv(dds_cons_res, "dds_cons.csv")

# ---- DESeq2 analysis: PI ----
mmrf_cln_per_pt$PI <- factor(mmrf_cln_per_pt$PI, levels = unique(mmrf_cln_per_pt$PI))
dds_tmp_PI <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts[, -1], colData = mmrf_cln_per_pt, design = ~PI)

# removing low-expressed genes
PI_exp_genes <- rowSums(counts(dds_tmp_PI)) >= 5
dds_PI <- dds_tmp_PI[PI_exp_genes,]

dds_PI <- DESeq(dds_PI)
dds_PI_res <- results(dds_PI, contrast = c("PI", "0", "1"), alpha = 0.05)
summary(dds_PI_res)

write.csv(dds_PPI_res, "dds_PI.csv")

# ---- DESeq2 analysis: IMID ----
mmrf_cln_per_pt$IMID <- factor(mmrf_cln_per_pt$IMID, levels = unique(mmrf_cln_per_pt$IMID))
dds_tmp_IMID <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts[, -1], colData = mmrf_cln_per_pt, design = ~IMID)

# removing low-expressed genes
IMID_exp_genes <- rowSums(counts(dds_tmp_IMID)) >= 5
dds_IMID <- dds_tmp_IMID[IMID_exp_genes,]

dds_IMID <- DESeq(dds_IMID)
dds_IMID_res <- results(dds_IMID, contrast = c("IMID", "0", "1"), alpha = 0.05)
summary(dds_IMID_res)

write.csv(dds_IMID_res, "dds_IMID.csv")

# ---- DESeq2 analysis: PI + IMID ----
mmrf_cln_per_pt$PI_IMID <- factor(mmrf_cln_per_pt$PI_IMID, levels = unique(mmrf_cln_per_pt$PI_IMID))
dds_tmp_PI_IMID <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts[, -1], colData = mmrf_cln_per_pt, design = ~PI_IMID)

# removing low-expressed genes
PI_IMID_exp_genes <- rowSums(counts(dds_tmp_PI_IMID)) >= 5
dds_PI_IMID <- dds_tmp_PI_IMID[PI_IMID_exp_genes,]

dds_PI_IMID <- DESeq(dds_PI_IMID)
dds_PI_IMID_res <- results(dds_PI_IMID, contrast = c("PI_IMID", "0", "1"), alpha = 0.05)
summary(dds_PI_IMID_res)

write.csv(dds_PI_IMID_res, "dds_PI_IMID.csv")

# ---- DESeq2 analysis: mAb ----
mmrf_cln_per_pt$mAb <- factor(mmrf_cln_per_pt$mAb, levels = unique(mmrf_cln_per_pt$mAb))
dds_tmp_mAb <- DESeqDataSetFromMatrix(countData = mmrf_cpm_counts[, -1], colData = mmrf_cln_per_pt, design = ~mAb)

# removing low-expressed genes
mAb_exp_genes <- rowSums(counts(dds_tmp_mAb)) >= 5
dds_mAb <- dds_tmp_mAb[mAb_exp_genes,]

dds_mAb <- DESeq(dds_mAb)
dds_mAb_res <- results(dds_mAb, contrast = c("mAb", "0", "1"), alpha = 0.05)
summary(dds_mAb_res)

write.csv(dds_mAb_res, "dds_mAb.csv")

