# !/usr/bin/r

## file: ratios.R
## last update: 29-07-2025


library(data.table)
library(ggplot2)
library(tidyverse)

source("C:/Users/Dell/Desktop/git_projects/TP53/scripts/fun/utils.R")


# ---- Main ----
# setting the env
wd <- setwd("C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")
outDir <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/visualization/"

# loading input
mmrf_tp53_per_pt <- fread("transcript_based/mmrf_tp53_per_pt.txt")


# ---- Data manipulation ----
# selecting isoforms exp
mmrf_exp_per_iso <- mmrf_tp53_per_pt %>%
  select(contains("p53")) %>%
  select(!contains("exp"))
mmrf_exp_per_iso <- mmrf_exp_per_iso + 1 #adding pseudocounts

# computing ratio
mmrf_tmp_ratio_per_iso <- mmrf_exp_per_iso %>%
  mutate(ratio_p53β_FL = p53β/p53_FL, 
         ratio_p53γ_FL = p53γ/p53_FL, 
         ratio_Δ40p53α_FL = Δ40p53α/p53_FL, 
         ratio_Δ133p53α_FL = Δ133p53α/p53_FL, 
         ratio_Δ133p53β_FL = Δ133p53β/p53_FL, 
         ratio_Δ133p53γ_FL = Δ133p53γ/p53_FL)

# computing cut-off value
mmrf_ratio <- as.data.frame(t(mmrf_ratio_per_iso %>%
  select(starts_with("ratio"))))

ratio_cut_offs <- defineTh(mmrf_ratio, "median")

# applying cut-off
mmrf_ratio_per_iso <- mmrf_tmp_ratio_per_iso %>%
  mutate(ratio_p53β_FL_exp = ifelse(ratio_p53β_FL >= ratio_cut_offs[["ratio_p53β_FL"]], "high_ratio", "low_ratio"), 
         ratio_p53γ_FL_exp = ifelse(ratio_p53γ_FL >= ratio_cut_offs[["ratio_p53γ_FL"]], "high_ratio", "low_ratio"),
         ratio_Δ40p53α_FL_exp = ifelse(ratio_Δ40p53α_FL >= ratio_cut_offs[["ratio_Δ40p53α_FL"]], "high_ratio", "low_ratio"),
         ratio_Δ133p53α_FL_exp = ifelse(ratio_Δ133p53α_FL >= ratio_cut_offs[["ratio_Δ133p53α_FL"]], "high_ratio", "low_ratio"),
         ratio_Δ133p53β_FL_exp = ifelse(ratio_Δ133p53β_FL >= ratio_cut_offs[["ratio_Δ133p53β_FL"]], "high_ratio", "low_ratio"),
         ratio_Δ133p53γ_FL_exp = ifelse(ratio_Δ133p53γ_FL >= ratio_cut_offs[["ratio_Δ133p53γ_FL"]], "high_ratio", "low_ratio"))

# binding
mmrf_lvl_rt_per_pt <- cbind(mmrf_tp53_per_pt %>% select(resp_sh, ends_with("exp")), mmrf_ratio_per_iso %>% select(ends_with("exp"))) %>%
  select(-p53_FL_exp) %>%
  filter(!resp_sh=="")
mmrf_lvl_rt_per_pt$resp_sh <- factor(mmrf_lvl_rt_per_pt$resp_sh, levels = c("PD", "SD", "PR", "VGPR", "CR", "sCR")) #c(""): should include those that do not resp?


# ---- Plotting ----
stackedBarRatio <- function(df, exp_col, ratio_col){
  return(df %>%
           select(resp_sh, !!sym(exp_col), !!sym(ratio_col)) %>%
           count(resp_sh, !!sym(exp_col), !!sym(ratio_col), name = "count") %>%
           ggplot(aes(x = !!sym(ratio_col), y = count, fill = !!sym(exp_col))) +
            geom_bar(position = "fill", stat = "identity") +
            facet_wrap(~resp_sh, nrow = 2) +
            labs(x = "", y = "") +
            scale_fill_manual(values = c("high" = "#B2182BFF", "low" = "#2166ACFF", "absent" = "darkgrey")) +
            theme_minimal()
         )
}


# p53β
mmrf_p53β_ratio <- stackedBarRatio(mmrf_lvl_rt_per_pt, "p53β_exp", "ratio_p53β_FL_exp") +
  labs(title = "p53β_exp/ratio_p53β_FL_exp")
ggsave(paste0(outDir, "p53β_exp_ratio.png"), mmrf_p53β_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# p53γ
mmrf_p53γ_ratio <- stackedBarRatio(mmrf_lvl_rt_per_pt, "p53γ_exp", "ratio_p53γ_FL_exp") +
  labs(title = "p53γ_exp/ratio_p53γ_FL_exp")
ggsave(paste0(outDir, "p53γ_exp_ratio.png"), mmrf_p53γ_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ40p53α
mmrf_Δ40p53α_ratio <- stackedBarRatio(mmrf_lvl_rt_per_pt, "Δ40p53α_exp", "ratio_Δ40p53α_FL_exp") +
  labs(title = "Δ40p53α_exp/ratio_Δ40p53α_FL_exp")
ggsave(paste0(outDir, "Δ40p53α_exp_ratio.png"), mmrf_Δ40p53α_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53α
mmrf_Δ133p53α_ratio <- stackedBarRatio(mmrf_lvl_rt_per_pt, "Δ133p53α_exp", "ratio_Δ133p53α_FL_exp") +
  labs(title = "Δ133p53α_exp/ratio_Δ133p53α_FL_exp")
ggsave(paste0(outDir, "Δ133p53α_exp_ratio.png"), mmrf_Δ133p53α_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53β
mmrf_Δ133p53β_ratio <- stackedBarRatio(mmrf_lvl_rt_per_pt, "Δ133p53β_exp", "ratio_Δ133p53β_FL_exp") +
  labs(title = "Δ133p53β_exp/ratio_Δ133p53β_FL_exp")
ggsave(paste0(outDir, "Δ133p53β_exp_ratio.png"), mmrf_Δ133p53β_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53γ
mmrf_Δ133p53γ_ratio <- stackedBarRatio(mmrf_lvl_rt_per_pt, "Δ133p53γ_exp", "ratio_Δ133p53γ_FL_exp") +
  labs(title = "Δ133p53γ_exp/ratio_Δ133p53γ_FL_exp")
ggsave(paste0(outDir, "Δ133p53γ_exp_ratio.png"), mmrf_Δ133p53γ_ratio, height = 10, width = 18, dpi = 400, bg = "white")

