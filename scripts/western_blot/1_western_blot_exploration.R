#!/usr/bin/r

## file: 1_western_blot_exploration
## aim: analyzed tp53 western-blot-derived protein info
## last update: 06-08-2026


# installing RNA-seq libraries
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install(c("drawProteins", "pRoloc", "pRolocGUI"))


# ---- Main ----
library(data.table)
library(DEqMS)
library(drawProteins)
library(ggplot2)
library(limma)
library(readxl)
library(tidyverse)

# custom function 
source("C:/Users/violameixian.vuong2/Desktop/git-projects/p53-isoforms/scripts/fun/utils.R")

# global inputs
wd <- "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/"
setwd(wd)

# WB data
pp_tp53 <- read_xlsx("original_db/Tabella pz TP53 (5).xlsx", sheet = "Densitometry") 
pts <- pp_tp53$MPC

# TP53 isoforms
tp53_df <- fread(paste0(wd, "tp53_isoforms.tsv"))

# ---- Pre-processing ----
wb_tmp_tp53 <- as.data.frame(pp_tp53) %>%
  select(`p53 FL (53 Kd)`, `p53β/γ (45 Kd)`, `Δ40p53α (42 Kd)`, `Δ40p53β (39 Kd)`, `Δ133p53α (35 Kd)`, `Δ160p53α (32 Kd)`, `Δ133p53 β (25 Kd)`) 
wb_tmp_tp53 <- as.data.frame(sapply(1:ncol(wb_tmp_tp53), function(x) as.numeric(wb_tmp_tp53[, x]))) %>% #convert to numeric values
  rename_with(~c("P04637_1", "P04637_2/3", "P04637_4", "P04637_5", "P04637_7", "A0A087X1Q1", "P04637_8")) 
wb_tp53 <- wb_tmp_tp53 %>% mutate(MPC = pts, .before = P04637_1) #adding pts

# ---- Normalization ----
## preparing for log normalization: adding pseudocount (+1)
wb_log_tp53 <- wb_tp53[, -1] + 1
wb_log_tp53 <- log2(wb_log_tp53) %>% na.omit() #MPC_2312 is NA
boxplot(wb_log_tp53, las=2) #suboptimalcheck for median centering


## limma-based data normalization
qt_tp53 <- normalizeBetweenArrays(wb_tp53[,-1])
plotDensities(qt_tp53) #check

# ---- DEA ----
## creating design matrix 
design <- model.matrix(~factor(c(1,2,3,4,5,6,7)))
colnames(design) <- c("P04637_1", "P04637_2", "P04637_4", "P04637_5", "P04637_7", "A0A087X1Q1", "P04637_8")

fit <- lmFit(qt_tp53, design)
fit2 <- contrasts.fit(fit, makeContrasts("P04637_1VSall" = P04637_1 - (P04637_2-P04637_4-P04637_5-P04637_7-A0A087X1Q1-P04637_8),
                                         #"P04637_2VSall" = P04637_2 - (P04637_1-P04637_4-P04637_5-P04637_7-A0A087X1Q1-P04637_8), 
                                         #"P04637_4VSall" = P04637_2 - (P04637_1-P04637_2-P04637_5-P04637_7-A0A087X1Q1-P04637_8), 
                                         #"P04637_5VSall" = P04637_4 - (P04637_1-P04637_2-P04637_4-P04637_7-A0A087X1Q1-P04637_8),
                                         #"P04637_7VSall" = P04637_7 - (P04637_1-P04637_2-P04637_4-P04637_5-A0A087X1Q1-P04637_8), 
                                         #"A0A087X1Q1VSall" = A0A087X1Q1 - (P04637_1-P04637_2-P04637_4-P04637_5-P04637_7-P04637_8),
                                         #"P04637_8VSall" = P04637_8 - (P04637_1-P04637_2-P04637_4-P04637_5-P04637_7-A0A087X1Q1), 
                                         levels = design))
fit2 <- eBayes(fit2)

qt_tp53 <- as.data.frame(qt_tp53) %>% mutate(MPC = pts)

test_diff(qt_tp53, type = "all")


# ---- Visualization ----
tp53_ann <- feature_to_dataframe(get_features("P04637"))
draw_domains(draw_canvas(tp53_ann), tp53_ann)


