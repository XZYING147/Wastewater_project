rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/FigureS4")

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

top_ten_df <- read.delim("rpkm.type_sample_composition.txt",sep="\t")
rownames(top_ten_df) <- top_ten_df[, 1]
top_ten_df <- top_ten_df[, -1]
top_ten_df <- as.data.frame(t(top_ten_df))
top_ten_df_with_row_names <- cbind(RowName = rownames(top_ten_df), top_ten_df)

long_df <- pivot_longer(top_ten_df_with_row_names, 
                        cols = -RowName, 
                        names_to = "genus", 
                        values_to = "value")

stack_df <- merge(long_df, samples_metadata, by.x = "RowName", by.y = "Sample_ID")
stack_df_overall <- stack_df %>%
  mutate(Sampling_year_month = "Overall")  # 将所有数据的月份标记为“Overall”
stack_df_extended <- bind_rows(stack_df, stack_df_overall)

stack_df_extended_averaged <- stack_df_extended %>%
  group_by(Sampling_type, genus,Sampling_year_month) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

p1_stack <- ggplot(data = stack_df_extended_averaged, aes(Sampling_type,value,fill=genus))+
  geom_bar(position="stack", stat = "identity",width=0.8, color = "white", size = 0.03)+ 
  facet_wrap(~ Sampling_year_month, nrow = 1, strip.position = "top") + 
  theme_linedraw() + 
  ylab("ARG RPKM") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15),  # x轴标签45度倾斜并末端对齐
    strip.background = element_rect(fill = "black"), # 分面标题背景设为黑色
    strip.text = element_text(color = "white", face = "bold", size = 20),  # 分面标题文字设为白色
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_x_discrete(labels = c("RES", "HOS", "WTP", "FLI")) +
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) +
  scale_fill_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                               "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                               "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                               "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) +
  xlab(element_blank()) 
p1_stack
ggsave("stack_resistome.pdf",p1_stack,device = "pdf",width = 12.5,height = 15)
