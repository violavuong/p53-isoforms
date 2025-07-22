#!/usr/bin/r

## file: visualization.R
## last update: 22-07-2025

# ---- Main ----
library(ComplexHeatmap)
library(data.table)
library(EnhancedVolcano)
library(ggplot2)
library(ggtranscript)
library(gridExtra)
library(patchwork)
library(readxl)
library(reshape2)
library(tidyverse)


wd <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/MM group - Vuong_Viola_Meixian/TP53/"
setwd(wd)

source("C:/Users/Dell/Desktop/git_projects/TP53/scripts/fun/utils.R")

# global input
mmrf_tp53_per_pt <- fread("mmrf_tp53_per_pt.txt")
mmrf_tp53_ann <- read_xlsx("mmrf_tp53.xlsx", sheet = "annotation")


# ---- cohort plot ----
demog_df <- mmrf_tp53_per_pt %>%
  mutate(age_group = ifelse(age > 65, "old", "young"), 
         gender = ifelse(gender==1, "male", "female")) %>%
  select(age_group, gender, ISS) #gender: 1 - male, 2 - female
demog_df$ISS <- as.character(demog_df$ISS)

demog_palette <- list(age_group = c(young = "#FED789FF", old = "#A4BED5FF"),
                      gender = c(male = "#6996E3FF", female = "#F28AAAFF"),
                      ISS = c("1" = "#96D2A4FF", "2" = "#6DBC90FF", "3" = "#266B6EFF"))

demog_ann <- HeatmapAnnotation(age = demog_df$age_group, gender = demog_df$gender, ISS = demog_df$ISS, 
                               col = demog_palette, which = "column", height = unit(4, "cm"),
                               #n_high = anno_barplot(num_high_isoforms_per_patient, border = FALSE,
                              #   gp = gpar(fill = "#8da0cb"),
                              #   axis_param = list(direction = "normal"),
                              #   height = unit(2.5, "cm") # Height for this bar plot
                               #),
                               annotation_name_side = "left")

tp53_expr_matrix <- as.matrix(t(mmrf_tp53_per_pt %>% select(ends_with("group"))))
tp53_color_fun <- c("low_exp" = "#4575b4", "high_exp" = "#d73027")

Heatmap(
  tp53_expr_matrix,
  name = "TP53 Expression", # Name for the legend
  col = tp53_color_fun,     # Apply the custom color scale
  
  # Clustering
  cluster_rows = TRUE,       
  cluster_columns = TRUE,    
  show_row_names = FALSE,    
  show_column_names = FALSE,  
  
  # Heatmap title
  column_title = "TP53 Isoform Expression across Patients",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  row_title = "isoform",
  row_title_gp = gpar(fontsize = 14),
  
  top_annotation = demog_ann, 
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "black", fill = fill, lwd = 0.5)) 
  },
  row_split = rownames(tp53_expr_matrix), 
  row_gap = unit(1, "mm"),
  border = FALSE
)


tp53_expr_matrix <- ifelse(tp53_expr_matrix=="low_exp", 0, 1)
pheat_col <- c("#4575b4", "#d73027")

pheatmap(tp53_expr_matrix,
         color = pheat_col, # Your custom colors for 0 and 1
         breaks = c(0, 1), # Ensure colors are mapped correctly to 0 and 1
         cellwidth = 12,
         cellheight = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,  # Show isoform names
         show_colnames = FALSE)

## heatmap with therapy and maybe ISS
mmrf_hmap_df <- as.data.frame(t(mmrf_l2_tp53[, -1])) %>%
  rownames_to_column("pt") %>%
  left_join(thr_df %>% select(PUBLIC_ID, D_PT_age, D_PT_gender, D_PT_iss, thername, condition), by = c("pt"="PUBLIC_ID")) %>%
  rename_with(~c(tp53_df$ensembl_ENST[-c(5,7)]), c(2,3,4,5,6,7))

hmap_data <- mmrf_hmap_df %>%
  select(pt, ENST00000269305, ENST00000420246, ENST00000455263, ENST00000504937, ENST00000510385, condition) %>%
  pivot_longer(cols = -c(pt, condition), names_to = "isoform", values_to = "value") %>%
  pivot_wider(names_from = isoform, values_from = value) %>%
  column_to_rownames("pt")

Heatmap(hmap_data, row_split = hmap_data$condition, cluster_rows = F, cluster_columns = F)


## testing oncoprint
mmrf_med_hmap_df <- as.data.frame(t(mmrf_tp53_median_class[, -755])) %>%
  rownames_to_column("pt") %>%
  left_join(thr_df %>% select(PUBLIC_ID, D_PT_age, D_PT_gender, D_PT_iss, thername, condition), by = c("pt"="PUBLIC_ID")) %>%
  rename_with(~c(tp53_df$ensembl_ENST[-c(4,5,7)]), c(2,3,4,5,6))

onco_data <- mmrf_med_hmap_df[, -755] %>%
  select(pt, ENST00000269305, ENST00000420246, ENST00000455263, ENST00000504937, ENST00000510385) %>%
  #pivot_longer(cols = -c(pt, condition), names_to = "isoform", values_to = "value") %>%
  #pivot_wider(names_from = isoform, values_from = value) %>%
  column_to_rownames("pt") %>%
  as.matrix()

hmap_ha <- rowAnnotation(age = mmrf_med_hmap_df[, -755]$D_PT_age, 
                         gender = as.factor(mmrf_med_hmap_df[, -755]$D_PT_gender), 
                         iss = as.factor(mmrf_med_hmap_df[, -755]$D_PT_iss), 
                         col = list(age = colorRamp2(c(65, max(mmrf_med_hmap_df[, -755]$D_PT_age, na.rm = T)), c("lightblue", "darkblue")),
                                    gender = c("1" = "#1f77b4", "2" = "#ff7f0e"),
                                    iss = c("1" = "#66c2a5", "2" = "#fc8d62", "3" = "#8da0cb")))

oncoPrint(onco_data, 
          get_type = function(x) {ifelse(x==1, "presence", ifelse(x==0, "absence", NA))}, 
          alter_fun = list(presence = function(x, y, w, h) {grid.rect(x, y, w, h, gp = gpar(fill = "#1f77b4", col = NA))}, 
                           absence = function(x, y, w, h) {grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey"))}), 
          col = colorRamp2(c(0, 1), c("white", "#1f77b4")),
          row_split = mmrf_med_hmap_df[, -755]$condition, # Split by condition
          row_title = "Condition",
          column_title = "Isoform Presence/Absence",
          heatmap_legend_param = list(title = "Presence/Absence", at = c(0, 1), labels = c("Absent", "Present")),
          top_annotation = HeatmapAnnotation(
            gender = factor(mmrf_med_hmap_df[, -755]$D_PT_gender),
            age = mmrf_med_hmap_df[, -755]$D_PT_age,
            iss = factor(mmrf_med_hmap_df[, -755]$D_PT_iss),
            col = list(gender = c("1" = "#1f77b4", "2" = "#ff7f0e"),
                       age = colorRamp2(c(65, max(mmrf_med_hmap_df[, -755]$D_PT_age)), c("lightblue", "darkblue")),
                       iss = c("1" = "#66c2a5", "2" = "#fc8d62", "3" = "#8da0cb")
            )
          ),
          show_column_names = TRUE,
          show_row_names = TRUE)


# ---- frequency stacked barplot ----
mmrf_tmp_tp53 <- mmrf_tp53_per_pt %>% select(p53_FL, p53β, p53γ, Δ40p53α, Δ133p53α, Δ133p53β, Δ133p53γ)
rownames(mmrf_tmp_tp53) <- mmrf_tp53_per_pt$PUBLIC_ID
mmrf_tp53 <- as.data.frame(t(mmrf_tmp_tp53))
colnames(mmrf_tp53) <- mmrf_tp53_per_pt$PUBLIC_ID

## adding a pseudo-count of 1
mmrf_tp53 <- mmrf_tp53 + 1

## defining cut-offs
medians <- defineTh(mmrf_tp53, "median") #median cut-off
perc_75 <- defineTh(mmrf_tp53, "perc_75") #75 percentile cut-off
perc_95 <- defineTh(mmrf_tp53, "perc_95") #95 percentile cut-off
cut_offs <- as.data.frame(cbind(medians, perc_75, perc_95))

## 1/0 conversion
mmrf_tp53_median_class <- binaryConversion(mmrf_tp53, cut_offs, th_col = 1, mmrf_tp53_per_pt$PUBLIC_ID, rownames(mmrf_tp53)) #by median value
mmrf_tp53_median_class$transcript <- factor(mmrf_tp53_median_class$transcript, levels = unique(mmrf_tp53_median_class$transcript))

mmrf_tp53_75perc_class <- binaryConversion(mmrf_tp53, cut_offs, th_col = 2, mmrf_tp53_per_pt$PUBLIC_ID, rownames(mmrf_tp53)) #by 75th percentile
mmrf_tp53_75perc_class$transcript <- factor(mmrf_tp53_75perc_class$transcript, levels = unique(mmrf_tp53_75perc_class$transcript))

mmrf_tp53_95perc_class <- binaryConversion(mmrf_tp53, cut_offs, th_col = 3, mmrf_tp53_per_pt$PUBLIC_ID, rownames(mmrf_tp53)) #by 95th percentile 
mmrf_tp53_95perc_class$transcript <- factor(mmrf_tp53_95perc_class$transcript, levels = unique(mmrf_tp53_95perc_class$transcript))


## median
mmrf_tp53_median_groups <- aggregateIsoforms(mmrf_tp53_median_class)
mmrf_tp53_median_groups$group <- factor(mmrf_tp53_median_groups$group, levels = unique(mmrf_tp53_median_groups$group))

### per isoform
freq_median <- freqBarplot(mmrf_tp53_median_class, "transcript") +
  labs(title = "Frequency per isoform - median cut-off", x = "isoform", y = "count", fill = "presence/absence")

### per group
group_median <- freqBarplot(mmrf_tp53_median_groups, "group") +
  labs(title = "Frequency per isoform group - median cut-off", x = "group", y = "count", fill = "presence/absence")

ggsave("mmrf_tp53_frequency_group_barplot.png", group_median, height = 10, width = 18, dpi = 400, bg = "white")


## 75th percentile
mmrf_tp53_75perc_groups <- aggregateIsoforms(mmrf_tp53_75perc_class)
mmrf_tp53_75perc_groups$group <- factor(mmrf_tp53_75perc_groups$group, levels = unique(mmrf_tp53_75perc_groups$group))

### per isoform
freq_75perc <- freqBarplot(mmrf_tp53_75perc_class, "transcript") +
  labs(title = "Frequency per isoform - 75th cut-off", x = "isoform", y = "count", fill = "presence/absence")

### per group
group_75perc <- freqBarplot(mmrf_tp53_75perc_groups, "group") +
  labs(title = "Frequency per isoform group - 75th cut-off", x = "group", y = "count", fill = "presence/absence")

ggsave("mmrf_tp53_frequency_barplot.png", 
       wrap_plots(freq_median, freq_75perc, group_median, group_75perc, nrow = 2) + plot_layout(guides = "collect"), 
       height = 10, width = 18, dpi = 400, bg = "white")


## 95th percentile
mmrf_tp53_95perc_groups <- aggregateIsoforms(mmrf_tp53_95perc_class)
mmrf_tp53_95perc_groups$group <- factor(mmrf_tp53_95perc_groups$group, levels = unique(mmrf_tp53_95perc_groups$group))

### per isoform
freq_95perc <- freqBarplot(mmrf_tp53_95perc_class, "transcript") +
  labs(title = "Frequency per isoform - 95th cut-off", x = "isoform", y = "count", fill = "presence/absence")

### per group
group_95perc <- freqBarplot(mmrf_tp53_95perc_groups, "group") +
  labs(title = "Frequency per isoform group - 95th cut-off", x = "group", y = "count", fill = "presence/absence")

# ---- transcripts/relative abundance/log relative abundance boxplot ----
## transcripts
mmrf_tp53_ann$transcript_name <- factor(mmrf_tp53_ann$transcript_name, levels = unique(mmrf_tp53_ann$transcript_name))

mmrf_tp53_ann_fig <- mmrf_tp53_ann %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
  geom_range(aes(fill = transcript_biotype)) +
  geom_intron(data = to_intron(mmrf_tp53_ann, "transcript_name"), 
              arrow.min.intron.length = 100000) +
  geom_text(data = add_exon_number(mmrf_tp53_ann, "transcript_name"), aes(x = (start + end) / 2, label = exon_number),
            size = 3.5, nudge_y = 0.4) +
  labs(title = "TP53 isoforms transcript", x = "Genomic coordinates", y = "Isoform") +
  scale_fill_manual(values = c("grey")) +
  scale_y_discrete(limits = rev) +
  theme() +
  theme_minimal()

ggsave("tp53_isoforms_transcripts.png", mmrf_tp53_ann_fig, height = 10, width = 18, dpi = 400, bg = "white")


## relative abundance
mmrf_tp53$transcript <- colnames(mmrf_tp53_per_pt)[2:7]
melt(mmrf_tp53) %>%
  ggplot(aes(x = transcript, y = value)) +
  geom_boxplot(aes(fill = transcript)) +
  coord_flip() +
  labs(title = "Boxplot per isoform", x = "isoform", y = "log2(CPM + 1)", fill = "isoform") +
  theme_minimal()


## log relative abundance
mmrf_log_tp53 <- log2(mmrf_tp53[, -755])
mmrf_l2_tp53 <- t(mmrf_l2_tp53[, -1])
colnames(mmrf_l2_tp53) <- transcripts

melt(mmrf_l2_tp53 %>% mutate(transcript = colnames(mmrf_tp53_per_pt)[2:7])) %>%
  ggplot(aes(x = Var2, y = value)) +
  geom_boxplot(aes(fill = Var2)) +
  coord_flip() +
  labs(title = "Boxplot per isoform", x = "isoform", y = "log2(CPM + 1)", fill = "isoform") +
  theme_minimal()




## gene transcripts genomic plots wrt tp53a
tp53_ann_no_canon <- tp53_ann_df  %>% filter(transcript_name != "TP53-201")
tp53_ann_rescaled <- to_diff(exons = tp53_ann_no_canon, ref_exons = tp53_ann_df, group_var = "transcript_name")

tp53_ann_df %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
  geom_range() +
  geom_intron(data = to_intron(tp53_ann_df, "transcript_name"), arrow.min.intron.length = 300) +
  geom_range(data = tp53_ann_rescaled, aes(fill = diff_type), alpha = 0.2)


# ---- volcano plots ----
# deseq inputs
dds_files <- list.files(path = wd, pattern = "^dds", recursive = TRUE, full.names = TRUE)
dds_dfs <- lapply(dds_files, function(x) fread(x, drop = "V1")) %>%
  setNames(str_remove_all(basename(dds_files), "dds_|\\.csv"))


findTopDEG <- function(dds_df){
  return(dds_df %>%
           filter(padj < 0.05) %>% 
           arrange(padj, desc(abs(log2FoldChange))) %>%
           select(transcript, log2FC = log2FoldChange, pvalue = padj) %>%
           mutate(log2FC = signif(log2FC, 3), 
                  pvalue = signif(pvalue, 3)) %>%
           head(n = 10))
}


createVolcano <- function(dds, dds_name, top_deg){
  return(EnhancedVolcano(dds, lab = dds$transcript, 
                         x = "log2FoldChange", y = "padj", 
                         title = toupper(paste0(dds_name)), subtitle = "",
                         selectLab = top_deg,
                         pCutoff = 0.05, FCcutoff = 1,
                         pointSize = 1.5, colAlpha = 0.7, labSize = 3.0,
                         cutoffLineType = 'dashed', cutoffLineCol = 'black'))
}




top_DEGs <- lapply(seq_along(dds_dfs), function(x) findTopDEG(dds_dfs[[x]]))
volcanos <- lapply(seq_along(dds_dfs), function(x) createVolcano(dds_dfs[[x]], names(dds_dfs), top_DEGs[[x]]$transcript))


volcanos[[1]] +
  annotation_custom(tableGrob(top_tbl, rows = NULL, theme = ttheme_minimal()), xmin = 4, ymin = 5)


