#!/usr/bin/r

## file: DRIMSeq.R
## last update: 04-06-2025

# installing required packages
BiocManager::install(c("DRIMSeq", "tximport"))

# ---- Main ----
library(BiocParallel)
library(data.table)
library(DRIMSeq)
library(GenomicFeatures)
library(tidyverse)
library(tximport)

setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/MM group - Vuong_Viola_Meixian/TP53/")

# ---- creating MMRF input file ----
pathDir <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - IA22/expression_estimates_transcript_based/"

transcripts <- fread(paste0(pathDir, "MMRF_CoMMpass_IA22_salmon_transcriptUnstrandedIgFiltered_counts.tsv"))$Transcript
mmrf_raw_counts <- filterMMRF(fread(paste0("MMRF_CoMMpass_IA22_salmon_transcriptUnstrandedIgFiltered_counts.tsv")), fread("mmrf_cln_per_pt.txt")) %>%
  mutate(transcript = transcripts)


# ---- creating tx2gene mapping file ----
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.114.gtf.gz")

# extracting mapping
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
colnames(tx2gene) <- c("TXNAME", "GENEID") #check that no version are present in ENSEMBL IDs


# ---- merge tx2gene mapping to MMRF data ----
mmrf_counts <- mmrf_raw_counts %>%
  left_join(tx2gene, by = c("transcript" = "TXNAME")) %>%
  dplyr::select(gene_id = GENEID, feature_id = transcript, everything())
mmrf_counts <- as.data.frame(na.omit(mmrf_counts))

write_tsv(mmrf_counts, "drim_mmrf_counts.txt")


# ---- custom function ----
runDRIMSeq <- function(counts, metadata, feature_col){
  dm_tmp <-  dmDSdata(counts = counts, samples = data.frame(sample_id = metadata$PUBLIC_ID, group = metadata$feature_col, stringsAsFactors = FALSE))
  
  # smallest number of samples to filter data
  min_num <- min(table(metadata$feature_col))
  
  dm <- dmFilter(dm_tmp,
                 min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                 min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)
  
  # running DRIMSeq
  register(SnowParam(workers = 4))
  model <- model.matrix(~group, data = DRIMSeq::samples(dm))
  dm <- dmPrecision(dm, design = model, BPPARAM = SnowParam(workers = 4))
  dm <- dmFit(dm, design = model, BPPARAM = SnowParam(workers = 4))
  dm <- dmTest(dm, coef = colnames(model)[2], BPPARAM = SnowParam(workers = 4))
  
  dm_res_gene <- results(dm, level = "gene")
  dm_res_gene$adj_pvalue <- p.adjust(dm_res_gene$pvalue, method = "BH")
  dm_res_transcript <- results(dm, level = "feature")
  dm_res_transcript$adj_pvalue <- p.adjust(dm_res_transcript$pvalue, method = "BH")
  
  # statistical significant genes/transcripts
  dm_genes <- dm_res_gene[dm_res_gene$adj_pvalue < 0.05, ]
  dm_transcripts <- dm_res_transcript[dm_res_transcript$adj_pvalue < 0.05, ]
  
  return(list(dm_genes, dm_transcripts))
}


# ---- ASCT ----
dm_tmp_asct <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$asct, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$asct))

dm_asct <- dmFilter(dm_tmp_asct,
                    min_samps_feature_expr = min_num, #transcripts must have at least 373 pt
                    min_feature_expr = 10, #min count for a feature in a pt
                    min_samps_feature_prop = min_num, #transcripts must have this proportion in at least 373 pt
                    min_feature_prop = 0.05, #min 5% proportion for a feature in a pt
                    min_gene_expr = 10, #gene must have at least 10 total counts across all pts
                    min_samps_gene_expr = 754 #gene must have min_gene_expr in ALL samples to be kept
                    )

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_asct))
dm_asct_precision <- dmPrecision(dm_asct, design = model, BPPARAM = SnowParam(workers = 4))
dm_asct_fit <- dmFit(dm_asct_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_asct <- dmTest(dm_asct_fit, coef = "group0", BPPARAM = SnowParam(workers = 4))

dm_asct_res_gene <- results(dm_asct, level = "gene")
dm_asct_res_gene$adj_pvalue <- p.adjust(dm_asct_res_gene$pvalue, method = "BH")
dm_asct_res_transcript <- results(dm_asct, level = "feature")
dm_asct_res_transcript$adj_pvalue <- p.adjust(dm_asct_res_transcript$pvalue, method = "BH")

# statistical significant genes/transcripts
dm_asct_genes <- dm_asct_res_gene[dm_asct_res_gene$adj_pvalue < 0.05, ]
dm_asct_transcripts <- dm_asct_res_transcript[dm_asct_res_transcript$adj_pvalue < 0.05, ]


# ---- n ASCT ----
dm_tmp_n_asct <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$n_asct, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$n_asct))

dm_n_asct <- dmFilter(dm_tmp_n_asct,
                      min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                      min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_n_asct))
dm_n_asct_precision <- dmPrecision(dm_n_asct, design = model, BPPARAM = SnowParam(workers = 4))
dm_n_asct_fit <- dmFit(dm_n_asct_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_n_asct <- dmTest(dm_n_asct_fit, coef = "group", BPPARAM = SnowParam(workers = 4))

dm_n_asct_res_gene <- results(dm_n_asct, level = "gene")
dm_n_asct_res_gene$adj_pvalue <- p.adjust(dm_n_asct_res_gene$pvalue, method = "BH")
dm_n_asct_res_transcript <- results(dm_n_asct, level = "feature")
dm_n_asct_res_transcript$adj_pvalue <- p.adjust(dm_n_asct_res_transcript$pvalue, method = "BH")

# statistical significant genes/transcripts
dm_n_asct_genes <- dm_n_asct_res_gene[dm_n_asct_res_gene$adj_pvalue < 0.05, ]
dm_n_asct_transcripts <- dm_n_asct_res_transcript[dm_n_asct_res_transcript$adj_pvalue < 0.05, ]


# ---- maintenance ----
dm_tmp_maint <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$maintenance, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$maintenance))

dm_maint <- dmFilter(dm_tmp_maint,
                     min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                     min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_maint))
dm_maint_precision <- dmPrecision(dm_maint, design = model, BPPARAM = SnowParam(workers = 4))
dm_maint_fit <- dmFit(dm_maint_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_maint <- dmTest(dm_maint_fit, coef = "group1", BPPARAM = SnowParam(workers = 4))

dm_maint_res_gene <- results(dm_maint, level = "gene")
dm_maint_res_gene$adj_pvalue <- p.adjust(dm_maint_res_gene$pvalue, method = "BH")
dm_maint_res_transcript <- results(dm_maint, level = "feature")
dm_maint_res_transcript$adj_pvalue <- p.adjust(dm_maint_res_transcript$pvalue, method = "BH")

# statistical significant genes/transcripts
dm_maint_genes <- dm_maint_res_gene[dm_maint_res_gene$adj_pvalue < 0.05, ]
dm_maint_transcripts <- dm_maint_res_transcript[dm_maint_res_transcript$adj_pvalue < 0.05, ]


# ---- maintenance lena ----
dm_tmp_maint_lena <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$maint_lena, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$maint_lena))

dm_maint_lena <- dmFilter(dm_tmp_maint_lena,
                          min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                          min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_maint_lena))
dm_maint_lena_precision <- dmPrecision(dm_maint_lena, design = model, BPPARAM = SnowParam(workers = 4))
dm_maint_lena_fit <- dmFit(dm_maint_lena_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_maint_lena <- dmTest(dm_maint_lena_fit, coef = "group", BPPARAM = SnowParam(workers = 4))

dm_maint_lena_res_gene <- results(dm_maint_lena, level = "gene")
dm_maint_lena_res_gene$adj_pvalue <- p.adjust(dm_maint_lena_res_gene$pvalue, method = "BH")
dm_maint_lena_res_transcript <- results(dm_maint_lena, level = "feature")
dm_maint_lena_res_transcript$adj_pvalue <- p.adjust(dm_maint_lena_res_transcript$pvalue, method = "BH")

# statistical significant genes/transcripts
dm_maint_lena_genes <- dm_maint_lena_res_gene[dm_maint_lena_res_gene$adj_pvalue < 0.05, ]
dm_maint_lena_transcripts <- dm_maint_lena_res_transcript[dm_maint_lena_res_transcript$adj_pvalue < 0.05, ]


# ---- consolidation ----
dm_tmp_cons <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$consolidation, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$consolidation))

dm_cons <- dmFilter(dm_tmp_cons,
                    min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                    min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_cons))
dm_cons_precision <- dmPrecision(dm_cons, design = model, BPPARAM = SnowParam(workers = 4))
dm_cons_fit <- dmFit(dm_cons_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_cons <- dmTest(dm_cons_fit, coef = "group1", BPPARAM = SnowParam(workers = 4))

dm_cons_res_gene <- results(dm_cons, level = "gene")
dm_cons_res_transcript <- results(dm_cons, level = "feature")


# ---- PI ----
dm_tmp_PI <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$PI, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$PI))

dm_PI <- dmFilter(dm_tmp_PI,
                  min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                  min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_PI))
dm_PI_precision <- dmPrecision(dm_PI, design = model, BPPARAM = SnowParam(workers = 4))
dm_PI_fit <- dmFit(dm_PI_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_PI <- dmTest(dm_PI_fit, coef = "group0", BPPARAM = SnowParam(workers = 4))

dm_PI_res_gene <- results(dm_PI, level = "gene")
dm_PI_res_transcript <- results(dm_PI, level = "feature")


# ---- IMID ----
dm_tmp_IMID <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$IMID, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$IMID))

dm_IMID <- dmFilter(dm_tmp_IMID,
                    min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                    min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_IMID))
dm_IMID_precision <- dmPrecision(dm_IMID, design = model, BPPARAM = SnowParam(workers = 4))
dm_IMID_fit <- dmFit(dm_IMID_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_IMID <- dmTest(dm_IMID_fit, coef = "group0", BPPARAM = SnowParam(workers = 4))

dm_IMID_res_gene <- results(dm_IMID, level = "gene")
dm_IMID_res_transcript <- results(dm_IMID, level = "feature")


# ---- PI_IMID ----
dm_tmp_PI_IMID <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$PI_IMID, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$PI_IMID))

dm_PI_IMID <- dmFilter(dm_tmp_PI_IMID,
                       min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                       min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_PI_IMID))
dm_PI_IMID_precision <- dmPrecision(dm_PI_IMID, design = model, BPPARAM = SnowParam(workers = 4))
dm_PI_IMID_fit <- dmFit(dm_PI_IMID_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_PI_IMID <- dmTest(dm_PI_IMID_fit, coef = "group0", BPPARAM = SnowParam(workers = 4))

dm_PI_IMID_res_gene <- results(dm_PI_IMID, level = "gene")
dm_PI_IMID_res_transcript <- results(dm_PI_IMID, level = "feature")


# ---- mAb ----
dm_tmp_mAb <- dmDSdata(counts = mmrf_counts, samples = data.frame(sample_id = mmrf_cln_per_pt$PUBLIC_ID, group = mmrf_cln_per_pt$mAb, stringsAsFactors = FALSE))

# smallest number of samples to filter data
min_num <- min(table(mmrf_cln_per_pt$mAb))

dm_mAb <- dmFilter(dm_tmp_mAb,
                   min_samps_feature_expr = min_num, min_feature_expr = 10, min_samps_feature_prop = min_num,
                   min_feature_prop = 0.05, min_gene_expr = 10, min_samps_gene_expr = 754)

# running DRIMSeq
register(SnowParam(workers = 4))
model <- model.matrix(~group, data = DRIMSeq::samples(dm_mAb))
dm_mAb_precision <- dmPrecision(dm_mAb, design = model, BPPARAM = SnowParam(workers = 4))
dm_mAb_fit <- dmFit(dm_mAb_precision, design = model, BPPARAM = SnowParam(workers = 4))
dm_mAb <- dmTest(dm_mAb_fit, coef = "group1", BPPARAM = SnowParam(workers = 4))

dm_mAb_res_gene <- results(dm_mAb, level = "gene")
dm_mAb_res_transcript <- results(dm_mAb, level = "feature")

tmp <- runDRIMSeq(mmrf_counts, mmrf_cln_per_pt, sym("mAb"))
