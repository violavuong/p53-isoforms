# !/usr/bin/r

## file: ratios.R
## last update: 24-07-2026

library(crosstable)
library(data.table)
library(flextable)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)
library(tidyverse)

source("C:/Users/violameixian.vuong2/Desktop/git-projects/p53-isoforms/scripts/fun/utils.R")


# ---- Main ----
# setting the env
wd <- setwd("C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/data/mmrf/")
outDir <- "C:/Users/Dell/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/visualization/"

# loading input
mmrf_tp53_per_pt <- fread("transcript_based/mmrf_tp53_per_pt.txt")


# ---- Data manipulation ----
# selecting isoforms exp
isoforms <- c("p53_FL$", "p53β$", "p53γ$", "Δ40p53α$", "Δ133p53α$", "Δ133p53β$", "Δ133p53γ$")
mmrf_exp_per_isoform <- mmrf_tp53_per_pt %>% select(PUBLIC_ID, matches(isoforms))

# computing ratio
mmrf_rt_per_isoform <- mmrf_exp_per_isoform %>%
  mutate(across(!c(PUBLIC_ID, p53_FL), .fns = ~ ./p53_FL, .names = "{.col}_FL_rt")) 


# ---- Beta vs FL ----
# 626 pts (40 without respsh)
mmrf_rt_beta_FL <- mmrf_rt_per_isoform %>%
  select(PUBLIC_ID, p53_FL, p53β, p53β_FL_rt) %>%
  filter(p53_FL!=0, p53β!=0)

beta_FL_th <- unname(defineTh(as.data.frame(t(mmrf_rt_beta_FL[, -1])), "median"))[3]

# applying cut-off
mmrf_rt_beta_FL_per_isoform <- mmrf_rt_beta_FL %>%
  mutate(rt_p53β_FL_exp = ifelse(p53β_FL_rt > beta_FL_th, "high_ratio", "low_ratio")) %>%
  left_join(mmrf_tp53_per_pt %>% select(PUBLIC_ID, resp_sh, p53_FL_exp, p53β_exp), by = "PUBLIC_ID")
mmrf_rt_beta_FL_per_isoform$resp_sh <- factor(mmrf_rt_beta_FL_per_isoform$resp_sh, levels = c("PD", "SD", "PR", "VGPR", "CR", "sCR"))

p53β_tbl <- as_flextable(crosstable(mmrf_rt_beta_FL_per_isoform %>% 
                                      filter(resp_sh!="") %>% 
                                      select(resp_sh, p53β_exp, rt_p53β_FL_exp), 
                                    by = c(rt_p53β_FL_exp, p53β_exp), label = FALSE, total = TRUE))

write.xlsx(mmrf_rt_beta_FL_per_isoform, "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/ratios_240726/mmrf_rt_beta_FL_per_pt.xlsx")


# ---- Gamma vs FL ----
# 181 pts (14 of which do not have resp_sh)
mmrf_rt_gamma_FL <- mmrf_rt_per_isoform %>%
  select(PUBLIC_ID, p53_FL, p53γ, p53γ_FL_rt) %>%
  filter(p53_FL!=0, p53γ!=0)

gamma_FL_th <- unname(defineTh(as.data.frame(t(mmrf_rt_gamma_FL[, -1])), "median"))[3]

# applying cut-off
mmrf_rt_gamma_FL_per_isoform <- mmrf_rt_gamma_FL %>%
  mutate(rt_p53γ_FL_exp = ifelse(p53γ_FL_rt > gamma_FL_th, "high_ratio", "low_ratio")) %>%
  left_join(mmrf_tp53_per_pt %>% select(PUBLIC_ID, resp_sh, p53_FL_exp, p53γ_exp), by = "PUBLIC_ID")
mmrf_rt_gamma_FL_per_isoform$resp_sh <- factor(mmrf_rt_gamma_FL_per_isoform$resp_sh, levels = c("PD", "SD", "PR", "VGPR", "CR", "sCR"))

p53γ_tbl <- as_flextable(crosstable(mmrf_rt_gamma_FL_per_isoform %>% 
                                      filter(resp_sh!="") %>% 
                                      select(resp_sh, p53γ_exp, rt_p53γ_FL_exp), 
                                    by = c(rt_p53γ_FL_exp, p53γ_exp), label = FALSE, total = TRUE))

write.xlsx(mmrf_rt_gamma_FL_per_isoform, "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/ratios_240726/mmrf_rt_gamma_FL_per_pt.xlsx")


# ---- Delta40 vs FL ----
# 9 pts (1 of which does not have resp_sh)
mmrf_rt_delta40_FL <- mmrf_rt_per_isoform %>%
  select(PUBLIC_ID, p53_FL, Δ40p53α, Δ40p53α_FL_rt) %>%
  filter(p53_FL!=0, Δ40p53α!=0)

delta40_FL_th <- unname(defineTh(as.data.frame(t(mmrf_rt_delta40_FL[, -1])), "median"))[3]

# applying cut-off
mmrf_rt_delta40_FL_per_isoform <- mmrf_rt_delta40_FL %>%
  mutate(rt_Δ40_FL_exp = ifelse(Δ40p53α_FL_rt > delta40_FL_th, "high_ratio", "low_ratio")) %>%
  left_join(mmrf_tp53_per_pt %>% select(PUBLIC_ID, resp_sh, p53_FL_exp, Δ40p53α_exp), by = "PUBLIC_ID")
mmrf_rt_delta40_FL_per_isoform$resp_sh <- factor(mmrf_rt_delta40_FL_per_isoform$resp_sh, levels = c("SD", "PR", "VGPR", "CR"))

p53Δ40_tbl <- as_flextable(crosstable(mmrf_rt_delta40_FL_per_isoform %>% 
                                      filter(resp_sh!="") %>% 
                                      select(resp_sh, Δ40p53α_exp, rt_Δ40_FL_exp), 
                                    by = c(rt_Δ40_FL_exp, Δ40p53α_exp), label = FALSE, total = TRUE))

write.xlsx(mmrf_rt_delta40_FL_per_isoform, "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/ratios_240726/mmrf_rt_delta40_FL_per_pt.xlsx")


# ---- Delta133a vs FL ----
# 509 pts (31 of which do not have resp_sh)
mmrf_rt_delta133a_FL <- mmrf_rt_per_isoform %>%
  select(PUBLIC_ID, p53_FL, Δ133p53α, Δ133p53α_FL_rt) %>%
  filter(p53_FL!=0, Δ133p53α!=0)

delta133a_FL_th <- unname(defineTh(as.data.frame(t(mmrf_rt_delta133a_FL[, -1])), "median"))[3]

# applying cut-off
mmrf_rt_delta133a_FL_per_isoform <- mmrf_rt_delta133a_FL %>%
  mutate(rt_Δ133a_FL_exp = ifelse(Δ133p53α_FL_rt > delta133a_FL_th, "high_ratio", "low_ratio")) %>%
  left_join(mmrf_tp53_per_pt %>% select(PUBLIC_ID, resp_sh, p53_FL_exp, Δ133p53α_exp), by = "PUBLIC_ID")
mmrf_rt_delta133a_FL_per_isoform$resp_sh <- factor(mmrf_rt_delta133a_FL_per_isoform$resp_sh, levels = c("PD", "SD", "PR", "VGPR", "CR", "sCR"))

p53Δ133a_tbl <- as_flextable(crosstable(mmrf_rt_delta133a_FL_per_isoform %>% 
                                        filter(resp_sh!="") %>% 
                                        select(resp_sh, Δ133p53α_exp, rt_Δ133a_FL_exp), 
                                      by = c(rt_Δ133a_FL_exp, Δ133p53α_exp), label = FALSE, total = TRUE))

write.xlsx(mmrf_rt_delta133a_FL_per_isoform, "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/ratios_240726/mmrf_rt_delta133a_FL_per_pt.xlsx")


# ---- Delta133b vs FL ----
# 11 pts
mmrf_rt_delta133b_FL <- mmrf_rt_per_isoform %>%
  select(PUBLIC_ID, p53_FL, Δ133p53β, Δ133p53β_FL_rt) %>%
  filter(p53_FL!=0, Δ133p53β!=0)

delta133b_FL_th <- unname(defineTh(as.data.frame(t(mmrf_rt_delta133b_FL[, -1])), "median"))[3]

# applying cut-off
mmrf_rt_delta133b_FL_per_isoform <- mmrf_rt_delta133b_FL %>%
  mutate(rt_Δ133b_FL_exp = ifelse(Δ133p53β_FL_rt > delta133b_FL_th, "high_ratio", "low_ratio")) %>%
  left_join(mmrf_tp53_per_pt %>% select(PUBLIC_ID, resp_sh, p53_FL_exp, Δ133p53β_exp), by = "PUBLIC_ID")
mmrf_rt_delta133b_FL_per_isoform$resp_sh <- factor(mmrf_rt_delta133b_FL_per_isoform$resp_sh, levels = c("SD", "PR", "VGPR", "CR"))

p53Δ133b_tbl <- as_flextable(crosstable(mmrf_rt_delta133b_FL_per_isoform %>% 
                                          filter(resp_sh!="") %>% 
                                          select(resp_sh, Δ133p53β_exp, rt_Δ133b_FL_exp), 
                                        by = c(rt_Δ133b_FL_exp, Δ133p53β_exp), label = FALSE, total = TRUE))

write.xlsx(mmrf_rt_delta133b_FL_per_isoform, "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/ratios_240726/mmrf_rt_delta133b_FL_per_pt.xlsx")


# ---- Saving all tbl ----
save_as_html("p53β" = p53β_tbl, "p53γ" = p53γ_tbl, "Δ40p53α" = p53Δ40_tbl, 
             "Δ133p53α" = p53Δ133a_tbl, "Δ133p53β" = p53Δ133b_tbl,
             path = "C:/Users/violameixian.vuong2/Alma Mater Studiorum Università di Bologna/PROJECT_TP53-isoforms - Documents/output/ratios_240726/exp_rt_frequency_tbl_240726.html")


# ---- Plotting: exp over ratio ----
stackedBarExpRatio <- function(df, exp_col, ratio_col){
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
mmrf_p53β_exp_rt <- stackedBarExpRatio(mmrf_rt_exp_per_isoform_per_resp, "p53β_exp", "rt_p53β_FL_exp") +
  labs(title = "p53β_exp/rt_p53β_FL_exp")
ggsave(paste0(outDir, "p53β_exp_ratio.png"), mmrf_p53β_exp_rt, height = 10, width = 18, dpi = 400, bg = "white")


# p53γ
mmrf_p53γ_exp_rt <- stackedBarExpRatio(mmrf_rt_exp_per_isoform_per_resp, "p53γ_exp", "rt_p53γ_FL_exp") +
  labs(title = "p53γ_exp/rt_p53γ_FL_exp")
ggsave(paste0(outDir, "p53γ_exp_ratio.png"), mmrf_p53γ_exp_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ40p53α
mmrf_Δ40p53α_exp_rt <- stackedBarExpRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ40p53α_exp", "rt_Δ40p53α_FL_exp") +
  labs(title = "Δ40p53α_exp/rt_Δ40p53α_FL_exp")
ggsave(paste0(outDir, "Δ40p53α_exp_ratio.png"), mmrf_Δ40p53α_exp_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53α
mmrf_Δ133p53α_exp_rt <- stackedBarExpRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53α_exp", "rt_Δ133p53α_FL_exp") +
  labs(title = "Δ133p53α_exp/rt_Δ133p53α_FL_exp")
ggsave(paste0(outDir, "Δ133p53α_exp_ratio.png"), mmrf_Δ133p53α_exp_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53β
mmrf_Δ133p53β_exp_rt <- stackedBarExpRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53β_exp", "rt_Δ133p53β_FL_exp") +
  labs(title = "Δ133p53β_exp/rt_Δ133p53β_FL_exp")
ggsave(paste0(outDir, "Δ133p53β_exp_ratio.png"), mmrf_Δ133p53β_exp_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53γ
mmrf_Δ133p53γ_exp_rt <- stackedBarExpRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53γ_exp", "rt_Δ133p53γ_FL_exp") +
  labs(title = "Δ133p53γ_exp/rt_Δ133p53γ_FL_exp")
ggsave(paste0(outDir, "Δ133p53γ_exp_ratio.png"), mmrf_Δ133p53γ_exp_rt, height = 10, width = 18, dpi = 400, bg = "white")


# ---- Plotting: exp-FL over ratio ----
stackedBarRatio <- function(df, exp_col, ratio_col){
  return(df %>%
           select(resp_sh, p53_FL_exp, !!sym(exp_col), !!sym(ratio_col)) %>%
           mutate(FL_isoform_exp = paste(p53_FL_exp, !!sym(exp_col), sep = "-")) %>%
           count(resp_sh, FL_isoform_exp, !!sym(ratio_col), name = "count") %>%
            ggplot(aes(x = !!sym(ratio_col), y = count, fill = FL_isoform_exp)) +
              geom_bar(position = "fill", stat = "identity") +
              facet_wrap(~resp_sh, nrow = 2) +
              labs(x = "", y = "") +
              scale_fill_brewer(palette = "YlGnBu") +
              theme_minimal()
  )
}

# p53β
mmrf_p53β_rt <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "p53β_exp", "rt_p53β_FL_exp") +
  labs(title = "p53FL_β_exp/rt_p53β_FL_exp")
ggsave(paste0(outDir, "p53β_rt.png"), mmrf_p53β_rt, height = 10, width = 18, dpi = 400, bg = "white")

# p53γ
mmrf_p53γ_rt <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "p53γ_exp", "rt_p53γ_FL_exp") +
  labs(title = "p53FL_γ_exp/rt_p53γ_FL_exp")
ggsave(paste0(outDir, "p53γ_rt.png"), mmrf_p53γ_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ40p53α
mmrf_Δ40p53α_rt <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ40p53α_exp", "rt_Δ40p53α_FL_exp") +
  labs(title = "Δ40p53α_FL_exp/rt_Δ40p53α_FL_exp")
ggsave(paste0(outDir, "Δ40p53α_rt.png"), mmrf_Δ40p53α_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53α
mmrf_Δ133p53α_rt <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53α_exp", "rt_Δ133p53α_FL_exp") +
  labs(title = "Δ133p53α_FL_exp/rt_Δ133p53α_FL_exp")
ggsave(paste0(outDir, "Δ133p53α_rt.png"), mmrf_Δ133p53α_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53β
mmrf_Δ133p53β_rt <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53β_exp", "rt_Δ133p53β_FL_exp") +
  labs(title = "Δ133p53β_FL_exp/rt_Δ133p53β_FL_exp")
ggsave(paste0(outDir, "Δ133p53β_rt.png"), mmrf_Δ133p53β_rt, height = 10, width = 18, dpi = 400, bg = "white")


# Δ133p53γ
mmrf_Δ133p53γ_rt <- stackedBarRatio(mmrf_rt_exp_per_isoform_per_resp, "Δ133p53γ_exp", "rt_Δ133p53γ_FL_exp") +
  labs(title = "Δ133p53γ_FL_exp/rt_Δ133p53γ_FL_exp")
ggsave(paste0(outDir, "Δ133p53γ_rt.png"), mmrf_Δ133p53γ_rt, height = 10, width = 18, dpi = 400, bg = "white")

