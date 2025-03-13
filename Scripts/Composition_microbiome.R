rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure1/Composition")

library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(stringr)

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

col_sums <- colSums(profile_bact_g_t_relative)
top_col_indices <- head(order(col_sums, decreasing = TRUE), 20)
top_col_names <- colnames(profile_bact_g_t_relative)[top_col_indices]
top_df <- profile_bact_g_t_relative[, top_col_names]
row_sums_top <- rowSums(top_df)
result_df <- data.frame(others = 1 - row_sums_top)
top_ten_df <- cbind(top_df, result_df)
top_ten_df_with_row_names <- cbind(RowName = rownames(top_ten_df), top_ten_df)

long_df <- pivot_longer(top_ten_df_with_row_names, 
                        cols = -RowName, 
                        names_to = "genus", 
                        values_to = "value")

composition_df <- merge(long_df, samples_metadata, by.x = "RowName", by.y = "Sample_ID")
composition_df_overall <- composition_df %>%
  mutate(Sampling_year_month = "Overall") 
composition_df_extended <- bind_rows(composition_df, composition_df_overall)

composition_df_extended_averaged <- composition_df_extended %>%
  group_by(Sampling_type, genus,Sampling_year_month) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

p1_composition <- ggplot(data = composition_df_extended_averaged, aes(Sampling_type,value*100,fill=genus))+
  geom_bar(stat="identity",position="stack",width=0.8, color = "white", size = 0.03)+ 
  facet_wrap(~ Sampling_year_month, nrow = 1, strip.position = "top") + 
  theme_linedraw() + 
  ylab("Relative abundance (%)") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15), 
    strip.background = element_rect(fill = "black"), 
    strip.text = element_text(color = "white", face = "bold", size = 20),
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_x_discrete(labels = c("RES", "HOS", "WTP", "FLI")) +
  scale_color_manual(values = c("#f4c0bd","#ead1d1","#d1eee9","#f7f4d7","#efd094","#bfdfd2","#e3edae","#ccdef4","#f6e2eb","#7ddacb","#f4ecb4","#c2e9e2","#f6dbf1","#a1c0c7","#d2d2d2","#f9e5b8","#e0dff1","#a9cbdf","#aec6da","#f5dfe1","#d6d7d8")) +
  scale_fill_manual(values = c("#f4c0bd","#ead1d1","#d1eee9","#f7f4d7","#efd094","#bfdfd2","#e3edae","#ccdef4","#f6e2eb","#7ddacb","#f4ecb4","#c2e9e2","#f6dbf1","#a1c0c7","#d2d2d2","#f9e5b8","#e0dff1","#a9cbdf","#aec6da","#f5dfe1","#d6d7d8")) +
  xlab(element_blank()) 
p1_composition
ggsave("compositon_microbiome.pdf",p1_composition,device = "pdf",width = 12.5,height = 15)