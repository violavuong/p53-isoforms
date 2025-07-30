# !/usr/bin/r

## file: ratios.R
## last update: 30-07-2025


library(data.table)
library(flextable)
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
isoforms <- c("p53_FL$", "p53β$", "p53γ$", "Δ40p53α$", "Δ133p53α$", "Δ133p53β$", "Δ133p53γ$")
mmrf_exp_per_isoform <- mmrf_tp53_per_pt %>% select(matches(isoforms))
mmrf_exp_per_isoform <- mmrf_exp_per_isoform + 1 #adding pseudocounts

# computing ratio
mmrf_rt_per_isoform <- mmrf_exp_per_isoform %>%
  mutate(across(!p53_FL, ~ ./p53_FL)) %>%
  select(-p53_FL) %>%
  rename_with(~c("rt_p53β_FL", "rt_p53γ_FL", "rt_Δ40p53α_FL", "rt_Δ133p53α_FL", "rt_Δ133p53β_FL", "rt_Δ133p53γ_FL"), c(1,2,3,4,5,6))

# computing cut-off value
rt_cut_offs <- unname(defineTh(as.data.frame(t(mmrf_rt_per_isoform)), "median"))

# applying cut-off
mmrf_tmp_rt_exp_per_isoform <- mmrf_rt_per_isoform %>%
  mutate(rt_p53β_FL_exp = ifelse(rt_p53β_FL > rt_cut_offs[1], "high_ratio", "low_ratio"), 
         rt_p53γ_FL_exp = ifelse(rt_p53γ_FL > rt_cut_offs[2], "high_ratio", "low_ratio"),
         rt_Δ40p53α_FL_exp = ifelse(rt_Δ40p53α_FL > rt_cut_offs[3], "high_ratio", "low_ratio"),
         rt_Δ133p53α_FL_exp = ifelse(rt_Δ133p53α_FL > rt_cut_offs[4], "high_ratio", "low_ratio"),
         rt_Δ133p53β_FL_exp = ifelse(rt_Δ133p53β_FL > rt_cut_offs[5], "high_ratio", "low_ratio"),
         rt_Δ133p53γ_FL_exp = ifelse(rt_Δ133p53γ_FL > rt_cut_offs[6], "high_ratio", "low_ratio"))

# binding
isoforms_exp <- c("p53β_exp", "p53γ_exp", "Δ40p53α_exp", "Δ133p53α_exp", "Δ133p53β_exp", "Δ133p53γ_exp")
mmrf_rt_exp_per_isoform <- cbind(mmrf_tp53_per_pt %>% select(resp_sh, matches(isoforms_exp)), 
                                 mmrf_tmp_rt_exp_per_isoform %>% select(ends_with("exp"))) 
mmrf_rt_exp_per_isoform[mmrf_rt_exp_per_isoform==""] <- NA

write_tsv(mmrf_rt_exp_per_isoform, "transcript_based/mmrf_rt_exp_per_isoform.txt")


# how many pts do not have a labelled response? 43 pts - analyses will be performed on 616 pts for which a trt response was added
mmrf_non_responders <- mmrf_rt_exp_per_isoform %>% filter(is.na(resp_sh))

mmrf_rt_exp_per_isoform_per_resp <- mmrf_rt_exp_per_isoform %>% filter(!is.na(resp_sh))
mmrf_rt_exp_per_isoform_per_resp$resp_sh <- factor(mmrf_rt_exp_per_isoform_per_resp$resp_sh, levels = c("PD", "SD", "PR", "VGPR", "CR", "sCR"))

write_tsv(mmrf_rt_exp_per_isoform_per_resp, "transcript_based/mmrf_rt_exp_per_isoform_per_resp.txt")


# ---- Frequency tbl ----
# p53β
p53β_df <- mmrf_rt_exp_per_isoform_per_resp %>% select(resp_sh, p53β_exp, rt_p53β_FL_exp)
p53β_tbl <- crosstable(p53β_df, by = c(rt_p53β_FL_exp, p53β_exp), label = FALSE, total = TRUE) %>%
  as_flextable(compact = TRUE, header_shown_n = 1:2)


# p53γ
p53γ_df <- mmrf_rt_exp_per_isoform_per_resp %>% select(resp_sh, p53γ_exp, rt_p53γ_FL_exp)
p53γ_tbl <- crosstable(p53γ_df, by = c(rt_p53γ_FL_exp, p53γ_exp), label = FALSE, total = TRUE) %>%
  as_flextable(compact = TRUE, header_shown_n = 1:2)


# Δ40p53α
Δ40p53α_df <- mmrf_rt_exp_per_isoform_per_resp %>% select(resp_sh, Δ40p53α_exp, rt_Δ40p53α_FL_exp)
Δ40p53α_tbl <- crosstable(Δ40p53α_df, by = c(rt_Δ40p53α_FL_exp, Δ40p53α_exp), label = FALSE, total = TRUE) %>%
  as_flextable(compact = TRUE, header_shown_n = 1:2)


# Δ133p53α
Δ133p53α_df <- mmrf_rt_exp_per_isoform_per_resp %>% select(resp_sh, Δ133p53α_exp, rt_Δ133p53α_FL_exp)
Δ133p53α_tbl <- crosstable(Δ133p53α_df, by = c(rt_Δ133p53α_FL_exp, Δ133p53α_exp), label = FALSE, total = TRUE) %>%
  as_flextable(compact = TRUE, header_shown_n = 1:2)


# Δ133p53β
Δ133p53β_df <- mmrf_rt_exp_per_isoform_per_resp %>% select(resp_sh, Δ133p53β_exp, rt_Δ133p53β_FL_exp)
Δ133p53β_tbl <- crosstable(Δ133p53β_df, by = c(rt_Δ133p53β_FL_exp, Δ133p53β_exp), label = FALSE, total = TRUE) %>%
  as_flextable(compact = TRUE, header_shown_n = 1:2)


# Δ133p53γ
Δ133p53γ_df <- mmrf_rt_exp_per_isoform_per_resp %>% select(resp_sh, Δ133p53γ_exp, rt_Δ133p53γ_FL_exp)
Δ133p53γ_tbl <- crosstable(Δ133p53γ_df, by = c(rt_Δ133p53γ_FL_exp, Δ133p53γ_exp), label = FALSE, total = TRUE) %>%
  as_flextable(compact = TRUE, header_shown_n = 1:2)


save_as_html("p53β" = p53β_tbl, "p53γ" = p53γ_tbl, "Δ40p53α" = Δ40p53α_tbl, 
             "Δ133p53α" = Δ133p53α_tbl, "Δ133p53β" = Δ133p53β_tbl, "Δ133p53γ" = Δ133p53γ_tbl, 
             path = "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/exp_rt_frequency_tbl.html")


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
mmrf_p53β_ratio <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "p53β_exp", "rt_p53β_FL_exp") +
  labs(title = "p53β_exp/rt_p53β_FL_exp")
ggsave(paste0(outDir, "p53β_exp_ratio.png"), mmrf_p53β_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# p53γ
mmrf_p53γ_ratio <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "p53γ_exp", "rt_p53γ_FL_exp") +
  labs(title = "p53γ_exp/rt_p53γ_FL_exp")
ggsave(paste0(outDir, "p53γ_exp_ratio.png"), mmrf_p53γ_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ40p53α
mmrf_Δ40p53α_ratio <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ40p53α_exp", "rt_Δ40p53α_FL_exp") +
  labs(title = "Δ40p53α_exp/rt_Δ40p53α_FL_exp")
ggsave(paste0(outDir, "Δ40p53α_exp_ratio.png"), mmrf_Δ40p53α_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53α
mmrf_Δ133p53α_ratio <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53α_exp", "rt_Δ133p53α_FL_exp") +
  labs(title = "Δ133p53α_exp/rt_Δ133p53α_FL_exp")
ggsave(paste0(outDir, "Δ133p53α_exp_ratio.png"), mmrf_Δ133p53α_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53β
mmrf_Δ133p53β_ratio <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53β_exp", "rt_Δ133p53β_FL_exp") +
  labs(title = "Δ133p53β_exp/rt_Δ133p53β_FL_exp")
ggsave(paste0(outDir, "Δ133p53β_exp_ratio.png"), mmrf_Δ133p53β_ratio, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53γ
mmrf_Δ133p53γ_ratio <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53γ_exp", "rt_Δ133p53γ_FL_exp") +
  labs(title = "Δ133p53γ_exp/rt_Δ133p53γ_FL_exp")
ggsave(paste0(outDir, "Δ133p53γ_exp_ratio.png"), mmrf_Δ133p53γ_ratio, height = 10, width = 18, dpi = 400, bg = "white")

