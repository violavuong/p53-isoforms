#!/usr/bin/r

# file: utils.R
# description: functions
# last update: 27-05-2025


aggregateIsoforms <- function(df){
  agg_df <- rbind(plyr::numcolwise(sum)(df[c(1,2,3),]), #TA
                  plyr::numcolwise(sum)(df[c(1,2,3,4),]), #long
                  plyr::numcolwise(sum)(df[c(5,6),]), #short
                  plyr::numcolwise(sum)(df[c(1,4,5),]), #alpha
                  plyr::numcolwise(sum)(df[c(2,6),]))
  agg_df[agg_df >= 1] <- 1 #if at least 1 then 1, otherwise 0
  agg_df <- agg_df %>% mutate(group = c("TA", "long", "short", "alpha", "beta"))
  
  return(agg_df)
}


binaryConversion <- function(df, th_df, th_col = c(1, 2), mmrf_pts, transcripts){
  df <- as.data.frame(t(sapply(1:nrow(df), function(x) ifelse(df[x, ] > th_df[x, th_col], 1, 0))))
  colnames(df) <- mmrf_pts
  df$transcript <- transcripts
  return(df)
}

createQQPlot <- function(df, transcript){
  return(df %>%
           ggplot(aes(sample = !!sym(transcript))) + 
            geom_qq() +
            geom_qq_line(color = "#476F84FF") +
            labs(title = paste0(transcript)) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5))
         )
}


defineTh <- function(df, th = c("median", "perc_75")){
  if (th=="median"){
    cut_off <- apply(df, 1, median, na.rm = T)
  } else if (th=="perc_75"){
    cut_off <- apply(df, 1, quantile, probs = 0.75)
  } else {
    message("wrong option.")
  }
  return(cut_off)
}


filterMMRF <- function(df, cln_df){
  df_per_pt <- df %>% select(transcript, grep("_1_BM_CD138pos", colnames(df)))
  colnames(df_per_pt) <- str_remove(colnames(df_per_pt), "_1_BM_CD138pos")
  df_per_pt <- df_per_pt %>% select(matches(cln_df$PUBLIC_ID))
  
  return(df_per_pt)
}


freqBarplot <- function(df, groups = c("transcript", "group")){
  return(pivot_longer(df, cols = -!!sym(groups), names_to = "pt", values_to = "value") %>%
           group_by(!!sym(groups), value) %>%
           summarise(count = n(), .groups = "drop") %>%
           ggplot(aes(x = !!sym(groups), y = count, fill = factor(value))) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = count), position = position_stack(0.5)) +
            scale_fill_manual(values = c("#B2182BFF", "#2166ACFF")) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                  axis.title = element_blank()) +
            coord_flip()
         )
}
