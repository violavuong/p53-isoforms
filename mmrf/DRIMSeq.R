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


# ---- DRIMSeq analysis ----
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
dm_asct_precision <- dmPrecision(dm_asct, 
                                 design = model.matrix(~group, data = DRIMSeq::samples(dm_asct)), 
                                 BPPARAM = BiocParallel::bpparam())
