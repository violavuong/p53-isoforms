#!/usr/bin/r

## file: exploration.R
## last update: 23-07-2025

# installing required packages
install.packages(c("ggridges", "gtsummary"))

# ---- Main ----
library(data.table)
library(ggh4x)
library(gmodels)
library(gtsummary)
library(patchwork)
library(readxl)
library(scales)
library(tidyverse)

wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")

# inputs
mmrf_cln_per_pt <- fread(paste0(wd, "/clinical_data/mmrf_cln_per_pt.txt"))
mmrf_cpm_per_pt <- fread(paste0(wd, "/transcript_based/mmrf_cpm_per_pt.txt"))
mmrf_tp53 <- read_excel("mmrf_tp53.xlsx", sheet = "isoforms")


# ---- filtering ----
## filtering for TP53 isoforms 
mmrf_tmp_tp53_cpm_per_pt <- mmrf_cpm_per_pt %>% 
  filter(transcript %in% mmrf_tp53$`ENSEMBL mRNA`)

## transposing the data
mmrf_tp53_cpm_per_pt <- as.data.frame(t(mmrf_tmp_tp53_cpm_per_pt %>% column_to_rownames(var = "transcript")))


# ---- exploration ----
## histogram with density curves
mmrf_tp53_cpm_per_pt %>%
  pivot_longer(cols = everything(), names_to = "transcript", values_to = "expression") %>%
  ggplot(aes(x = expression)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = NULL, fill = "#A4BED5FF", color = "black", alpha = 0.7) +
    geom_density(color = "#476F84FF", linewidth = 1) +
    facet_wrap(~ transcript, scales = "free") +
    labs(title = "distribution of TP53 transcripts expression", x = "expression level", y = "density/frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), 
          strip.text = element_text(face = "bold"))

ggsave("mmrf_tp53_cpm_histogram.png", p, height = 10, width = 15, dpi = 400, bg = "white")

## boxplots
mmrf_tp53_cpm_per_pt %>%
  pivot_longer(cols = everything(), names_to = "transcript", values_to = "expression") %>%
  ggplot(aes(x = transcript, y = expression, fill = transcript)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = "boxplots of TP53 transcripts expression", x = "transcript", y = "expression") +
    scale_fill_manual(values = c("#FED789FF", "#023743FF", "#72874EFF", "#476F84FF", "#A4BED5FF", "#453947FF")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))

ggsave("mmrf_tp53_cpm_boxplot.png", p, height = 8, width = 15, dpi = 400, bg = "white")


## Q-Q plots
qq_plots <- lapply(colnames(mmrf_tp53_cpm_per_pt), function(x) createQQPlot(mmrf_tp53_cpm_per_pt, x))
wrap_plots(qq_plots) +
  plot_annotation(title = "QQ-plots for normal distribution")

ggsave("mmrf_tp53_cpm_qq.png", p, height = 10, width = 15, dpi = 400, bg = "white")


## descriptive exploration
sink("mmrf_tp53_cpm_per_pt_summary.txt")
print(summary(mmrf_tp53_cpm_per_pt))
sink() 


# ---- categorization ----
median_per_tp53 <- mmrf_tp53_cpm_per_pt %>% summarise(across(everything(), ~median(.[. > 0])))

mmrf_tmp_tp53_per_pt <- mmrf_tp53_cpm_per_pt %>%
  rownames_to_column("PUBLIC_ID") %>%
  left_join(mmrf_cln_per_pt, by = "PUBLIC_ID") %>% # Join with your clinical data
  mutate(across(.cols = starts_with("ENST"), 
                .fns = ~ factor(ifelse(.x==0, "absent", ifelse(.x > median_per_tp53[[cur_column()]], "high", "low")), levels = c("absent", "low", "high")),
                .names = "{.col}_exp"), 
         resp_group = ifelse(respsh=="", NA, ifelse(respsh %in% c("sCR", "CR", "VGPR"), "responder", "non_responder")), .after = respsh)

mmrf_tp53_per_pt <- mmrf_tmp_tp53_per_pt %>%
  select(PUBLIC_ID, age, gender, ISS, ECOG, 
         ALB, Ca, creatinine, HE, LDH, LDH_level, serum_M_prot, urine_24h_M_prot, serum_B2_mglob, serum_PCR, light_chain, light_chain_type, PLT, CD138_perc,
         ther_end = therendy_max, ther_name = thername, ther_cat = thercat, resp, resp_sh = respsh, resp_group, 
         PI, IMID, PI_IMID, mAb, asct, n_asct, maintenance, maint_lena, consolidation,
         start_line_1, end_line_1, best_resp_dy_line_1, PFS_date, PFS_event, PFS_time, OS_date, OS_event, OS_time,
         p53_FL = ENST00000269305, p53_FL_exp = ENST00000269305_exp, 
         p53β = ENST00000420246, p53β_exp = ENST00000420246_exp, 
         p53γ = ENST00000455263, p53γ_exp = ENST00000455263_exp, 
         Δ40p53α = ENST00000610292, Δ40p53α_exp = ENST00000610292_exp, 
         Δ133p53α = ENST00000504937, Δ133p53α_exp = ENST00000504937_exp, 
         Δ133p53β = ENST00000510385, Δ133p53β_exp = ENST00000510385_exp, 
         Δ133p53γ = ENST00000504290, Δ133p53γ_exp = ENST00000504290_exp)
mmrf_tp53_per_pt$light_chain_type <- tolower(mmrf_tp53_per_pt$light_chain_type)
mmrf_tp53_per_pt$resp_group <- factor(mmrf_tp53_per_pt$resp_group, levels = unique(mmrf_tp53_per_pt$resp_group))

write_tsv(mmrf_tp53_per_pt, "transcript_based/mmrf_tp53_per_pt.txt")


# ---- Fisher/chi-square significance ----
isoforms <- mmrf_tp53$`ENSEMBL mRNA`

## fisher test: high/low-expressing pts vs responder/non-responders
### using crosstable
fisher_summary <- tibble(isoform = character(), fisher_p_value = numeric(), odds_ratio = numeric(),
                         CI_lower = numeric(), CI_upper = numeric(),
                         n_analyzed = integer(), n_tot = integer())

for (isoform in isoforms) {
  cat("Analyzing", isoform, "expression vs therapy response\n")
  
  cross_tbl <- CrossTable(x = mmrf_tp53_per_pt[[paste0(isoform, "_exp")]],
                          y = mmrf_tp53_per_pt$resp_group,
                          prop.chisq = FALSE, prop.t = FALSE, prop.r = TRUE, prop.c = TRUE, 
                          chisq = FALSE, fisher = TRUE,
                          dnn = c(isoform, "response"))
  
  # extract results, ct$fisher.ts is for two-sided test
  p_value <- odds_ratio <- CI_lower <- CI_upper <- NA
  
  if (!is.null(cross_tbl$fisher.ts)) {
    p_value <- cross_tbl$fisher.ts$p.value
    odds_ratio <- cross_tbl$fisher.ts$estimate
    CI_lower <- cross_tbl$fisher.ts$conf.int[1]
    CI_upper <- cross_tbl$fisher.ts$conf.int[2]
  } else if (!is.null(cross_tbl$fisher.lt)) {
    p_value <- cross_tbl$fisher.lt$p.value
    odds_ratio <- cross_tbl$fisher.lt$estimate
    CI_lower <- cross_tbl$fisher.lt$conf.int[1]
    CI_upper <- cross_tbl$fisher.lt$conf.int[2]
  } else if (!is.null(cross_tbl$fisher.gt)) {
    p_value <- cross_tbl$fisher.gt$p.value
    odds_ratio <- cross_tbl$fisher.gt$estimate
    CI_lower <- cross_tbl$fisher.gt$conf.int[1]
    CI_upper <- cross_tbl$fisher.gt$conf.int[2]
  }
  
  n_analyzed <- sum(cross_tbl$t)
  fisher_summary <- fisher_summary %>%
    add_row(isoform = isoform, fisher_p_value = p_value, odds_ratio = as.numeric(odds_ratio), 
            CI_lower = CI_lower, CI_upper = CI_upper, n_analyzed = n_analyzed, n_tot = 754)
}

write_tsv(fisher_summary, "fisher_summary.txt")

### using tbl_summary
mmrf_tp53_per_pt %>%
  select(ends_with("exp"), resp_group) %>%
  tbl_summary(by = resp_group,
              include = c(ends_with("exp")),
              missing = "ifany",
              missing_text = "(Unknown)") %>%
  add_p(test = list(all_categorical() ~ "fisher.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  add_q() %>%
  modify_header(stat_by = "**{level}**<br>N = {n}") %>%
  modify_caption(caption = "**Total number of patients = {nrow(mmrf_tp53_per_pt)}. \n NAs are excluded from statistical tests**")


## chi-square test
### using crosstable
chi_summary <- tibble(isoform = character(), chi_p_value = numeric(), stat = numeric(), dof = numeric(),
                      n_analyzed = integer(), n_tot = integer())

for (isoform in isoforms) {
  cat("Analyzing", isoform, "expression vs therapy response\n")
  
  cross_tbl <- CrossTable(x = mmrf_tp53_per_pt[[paste0(isoform, "_exp")]],
                          y = mmrf_tp53_per_pt$resp_group,
                          prop.chisq = FALSE, prop.t = FALSE, prop.r = TRUE, prop.c = TRUE, 
                          chisq = TRUE, fisher = FALSE,
                          expected = FALSE, 
                          dnn = c(isoform, "response"))
  
  # extract results, ct$fisher.ts is for two-sided test
  p_value <- stat <- dof <- NA
  
  if (!is.null(cross_tbl$chisq)) {
    p_value <- cross_tbl$chisq$p.value
    stat <- cross_tbl$chisq$statistic
    dof <- cross_tbl$chisq$parameter
  }
  
  n_analyzed <- sum(cross_tbl$t)
  chi_summary <- chi_summary %>%
    add_row(isoform = isoform, chi_p_value = p_value, stat = stat, dof = dof, n_analyzed = n_analyzed, n_tot = 754)
}

write_tsv(chi_summary, "chi_summary.txt")

### using tbl_summary
mmrf_tp53_per_pt %>%
  select(ends_with("exp"), resp_group) %>%
  tbl_summary(by = resp_group,
              include = c(ends_with("exp")),
              missing = "ifany",
              missing_text = "(Unknown)") %>%
  add_p(test = list(all_categorical() ~ "chisq.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  add_q() %>%
  modify_header(stat_by = "**{level}**<br>N = {n}") %>%
  modify_caption(caption = "**Total number of patients = {nrow(mmrf_tp53_per_pt)}. \n NAs are excluded from statistical tests**")


# ---- statistical significance visualization ----
## stacked barplot
fisher_summary <- fisher_summary %>%
  mutate(label = paste0(isoform, " Fisher p-value: ", format.pval(fisher_p_value))) 

mmrf_tp53_long <- mmrf_tp53_per_pt %>%
  filter(!is.na(resp_group)) %>%
  select(ends_with("exp"), resp_group) %>%
  pivot_longer(cols = ends_with("exp"), names_to = "isoform", values_to = "exp") %>%
  mutate(isoform = gsub("_exp", "", isoform)) %>%
  left_join(fisher_summary %>% select(isoform, label), by = "isoform") %>%
  mutate(label = factor(label, levels = fisher_summary$label), 
         isoform = factor(isoform, levels = unique(isoform)))

mmrf_tp53_long %>%
  ggplot(aes(x = resp_group, fill = exp)) +
    geom_bar(position = "fill", color = "black", width = 0.7) +
    facet_wrap(~label, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = c("#2166ACFF", "#B2182BFF")) +
    scale_y_continuous(labels = percent) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    coord_flip()

## association
par(mfrow = c(2, 3))

for(tp53_isoform in isoforms){
  isoform_df <- mmrf_tp53_long %>%
    filter(isoform==tp53_isoform)
  isoform_tbl <- table(isoform_df$resp_group, isoform_df$exp)
  assocplot(tbl, 
            main = unique(isoform_df$label), xlab = "Clinical response group", ylab = "Expression category", col = c("#2166ACFF", "#B2182BFF"))
}
