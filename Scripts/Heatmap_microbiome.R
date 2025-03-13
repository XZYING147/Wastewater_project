rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/FigureS2")

library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

samples_metadata <- read.csv("Sample_metadata.csv", fileEncoding = "GBK") %>%
  mutate(
    Sampling_time = as.Date(Sampling_time, format = "%Y/%m/%d"),
    Sampling_year_month = format(Sampling_time, "%Y-%m")
  )

profile_bact <- read.delim("metaphlan4_taxonomy.tsv")

profile_bact$min_level <- str_extract(profile_bact$ID, "(k__|p__|c__|o__|f__|g__|s__|t__)[^|]*$")

profile_bact_g <- subset(profile_bact, grepl("^g__", profile_bact$min_level))[, -ncol(profile_bact)]

transpose_and_set_rowname <- function(df) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
  df <- as.data.frame(t(df))
  return(df)
}

profile_bact_g_t <- transpose_and_set_rowname(profile_bact_g)

profile_bact_g_t_relative <- profile_bact_g_t / rowSums(profile_bact_g_t)

col_sums_result <- colSums(profile_bact_g_t_relative)
profile_bact_g_t_relative <- profile_bact_g_t_relative[, col_sums_result > 0]
profile_bact_g_t_relative <- profile_bact_g_t_relative[, colSums(profile_bact_g_t_relative) > 0]

profile_bact_g_t_relative$Sample <- rownames(profile_bact_g_t_relative)
heatmap_metadata_df <- merge(profile_bact_g_t_relative, samples_metadata, by.x = "Sample", by.y = "Sample_ID")
heatmap_metadata_df_sorted <- heatmap_metadata_df %>% arrange(Sampling_type)
col_sums <- colSums(heatmap_metadata_df_sorted[1:252,2:769])

top_100_index <- order(col_sums, decreasing = TRUE)[1:100]
top_100_matrix <- as.matrix(heatmap_metadata_df_sorted[,top_100_index+1])
rownames(top_100_matrix) <- heatmap_metadata_df_sorted$Sample
top100_df=data.frame(top_100_matrix)
top100_df$Sample <- rownames(top100_df)

samples_metadata_sorted <- samples_metadata %>% arrange(Sampling_type)
top100_df$Sample <- samples_metadata_sorted$Sampling_type

group_fli_rows <- top100_df %>% filter(Sample == "Flight")
group_hos_rows <- top100_df %>% filter(Sample == "HOS")
group_res_rows <- top100_df %>% filter(Sample == "RES")
group_wtp_rows <- top100_df %>% filter(Sample == "WTP")

na_rows_for_group_a <- tibble(Sample = rep(NA, 2))
na_rows_for_group_b <- tibble(Sample = rep(NA, 2))
na_rows_for_group_c <- tibble(Sample = rep(NA, 2))

heatmap_na_df <- bind_rows(group_fli_rows,na_rows_for_group_a,group_hos_rows,na_rows_for_group_b,group_res_rows,na_rows_for_group_c,group_wtp_rows)
colnames(heatmap_na_df) <- gsub("_unclassified", "", colnames(heatmap_na_df))
heatmap_na_df_renamed <- heatmap_na_df %>%
  rename_with(~ str_extract(., "(?<=_)[^_]+$"), contains("_"))

heatmap_na_df_renamed <- heatmap_na_df_renamed[, -ncol(heatmap_na_df_renamed)]
profile_bact_g_t_relative_heatmap <- log10(heatmap_na_df_renamed+1)
profile_bact_g_t_relative_heatmap_turn <- t(profile_bact_g_t_relative_heatmap)

group_fli_rows_2 <- samples_metadata_sorted %>% filter(Sampling_type == "Flight")
group_hos_rows_2 <- samples_metadata_sorted %>% filter(Sampling_type == "HOS")
group_res_rows_2 <- samples_metadata_sorted %>% filter(Sampling_type == "RES")
group_wtp_rows_2 <- samples_metadata_sorted %>% filter(Sampling_type == "WTP")

na_rows_for_group_a_2 <- tibble(Sampling_type = rep(NA, 2))
na_rows_for_group_b_2 <- tibble(Sampling_type = rep(NA, 2))
na_rows_for_group_c_2 <- tibble(Sampling_type = rep(NA, 2))

metadata_na_df <- bind_rows(group_fli_rows_2,na_rows_for_group_a_2,group_hos_rows_2,na_rows_for_group_b_2,group_res_rows_2,na_rows_for_group_c_2,group_wtp_rows_2)

annotation_cols<- data.frame(Group = metadata_na_df$Sampling_type,
                             Batch = metadata_na_df$Batch,
                             Month = metadata_na_df$Sampling_year_month,
                             Area = metadata_na_df$Sampling_area)

rownames(annotation_cols) <- colnames(profile_bact_g_t_relative_heatmap_turn)

Group_color <- c("#f0a8c2", "#71c7ea", "#aed2d2", "#ecd8aa","white") 
names(Group_color) <- c("RES", "HOS", "WTP", "Flight","NA")

Batch_color <- c("#9998ff","#f4a99b","#bd9aad","white")
names(Batch_color) <- c("B1","B2","B3","NA")

Month_color <- c("#94baf5","#d6c4e0","#e2e2b6","white") 
names(Month_color) <- c("2023-02","2023-03","2023-04","NA")

Area_color <- c("#f4c0bd","#ead1d1","#d1eee9","#f7f4d7","#efd094","#aec6da","#e3edae","white")
names(Area_color) <- c("HB","HL","JM","TA","XA","HC","SM","NA") 

ann_colors <- list(Group=Group_color,Batch=Batch_color,Month=Month_color,Area=Area_color)

p1_heatmap <- pheatmap(profile_bact_g_t_relative_heatmap_turn,
                       # scale="row",#对行进行归一化
                       color = colorRampPalette(c("#ebeae6","#b23c49" ))(30), 
                       annotation_col = annotation_cols,
                       annotation_colors = ann_colors,
                       fontsize_col = 1, 
                       cluster_rows = F,
                       cluster_cols = F,
                       show_rownames =T, 
                       show_colnames = F,
                       fontsize = 2.5,
                       cellwidth=2.5,
                       cellheight=2.5,
                       main = "Microbiome Genera Level")

ggsave("heatmap_microbiome.pdf",p1_heatmap,device = "pdf",width = 15,height = 20)
