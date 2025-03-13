rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure2/Stream")

library(tidyverse)
library(patchwork)

data <- read.table("arg_type_with_others.csv",header = T,sep = ",",row.names = 1) 

data_fli <- data[1:22,1:22]
data_hos <- data[1:22,23:38]
data_res <- data[1:22,39:93]
data_wtp <- data[1:22,94:252]

first_col_fli <- data_fli[, 1]
combined_col_fli <- c(first_col_fli, unlist(data_fli[, -1]))
result_df_fli <- as.data.frame(combined_col_fli)
result_df_fli$x <- rep(c(1:22), each = 22)
result_df_fli$group <- rep(c("Aminoglycoside","Antibacterial_fatty_acid","Bacitracin","Beta_lactam","Bleomycin","Chloramphenicol","Florfenicol","Fosfomycin","MLS","Multidrug","Mupirocin","Novobiocin","Pleuromutilin_tiamulin","Polymyxin","Quinolone","Rifamycin","Streptothricin","Sulfonamide","Tetracycline","Trimethoprim","Vancomycin","Others"), 22)

first_col_hos <- data_hos[, 1]
combined_col_hos <- c(first_col_hos, unlist(data_hos[, -1]))
result_df_hos <- as.data.frame(combined_col_hos)
result_df_hos$x <- rep(c(1:16), each = 22)
result_df_hos$group <- rep(c("Aminoglycoside","Antibacterial_fatty_acid","Bacitracin","Beta_lactam","Bleomycin","Chloramphenicol","Florfenicol","Fosfomycin","MLS","Multidrug","Mupirocin","Novobiocin","Pleuromutilin_tiamulin","Polymyxin","Quinolone","Rifamycin","Streptothricin","Sulfonamide","Tetracycline","Trimethoprim","Vancomycin","Others"), 16)

first_col_res <- data_res[, 1]
combined_col_res <- c(first_col_res, unlist(data_res[, -1]))
result_df_res <- as.data.frame(combined_col_res)
result_df_res$x <- rep(c(1:55), each = 22)
result_df_res$group <- rep(c("Aminoglycoside","Antibacterial_fatty_acid","Bacitracin","Beta_lactam","Bleomycin","Chloramphenicol","Florfenicol","Fosfomycin","MLS","Multidrug","Mupirocin","Novobiocin","Pleuromutilin_tiamulin","Polymyxin","Quinolone","Rifamycin","Streptothricin","Sulfonamide","Tetracycline","Trimethoprim","Vancomycin","Others"), 55)

first_col_wtp <- data_wtp[, 1]
combined_col_wtp <- c(first_col_wtp, unlist(data_wtp[, -1]))
result_df_wtp <- as.data.frame(combined_col_wtp)
result_df_wtp$x <- rep(c(1:159), each = 22)
result_df_wtp$group <- rep(c("Aminoglycoside","Antibacterial_fatty_acid","Bacitracin","Beta_lactam","Bleomycin","Chloramphenicol","Florfenicol","Fosfomycin","MLS","Multidrug","Mupirocin","Novobiocin","Pleuromutilin_tiamulin","Polymyxin","Quinolone","Rifamycin","Streptothricin","Sulfonamide","Tetracycline","Trimethoprim","Vancomycin","Others"), 159)

p_fli <- ggplot(result_df_fli) +
  geom_area(aes(x, combined_col_fli, fill = group),
            position = "fill") +
  scale_fill_manual(values = rev(c("#f4c0bd", "#ead1d1", "#d1eee9", "#f7f4d7", "#efd094",
                                   "#bfdfd2", "#e3edae", "#ccdef4", "#f6e2eb", "#7fd8cc",
                                   "#f4ecb4", "#c2e9e2", "#d6d7d8", "#f6dbf1", "#a1c0c7",
                                   "#e8ddaf", "#d2d2d2", "#718991", "#b8b3d6", "#a9cbdf",
                                   "#f5dfe1", "#aec6da")))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("")+
  ylab("")+
  ggtitle("Flight")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        # axis.text.y = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), 
        legend.position = "none")
p_fli
p_hos <- ggplot(result_df_hos) +
  geom_area(aes(x, combined_col_hos, fill = group),
            position = "fill") +
  scale_fill_manual(values = rev(c("#f4c0bd", "#ead1d1", "#d1eee9", "#f7f4d7", "#efd094",
                                   "#bfdfd2", "#e3edae", "#ccdef4", "#f6e2eb", "#7fd8cc",
                                   "#f4ecb4", "#c2e9e2", "#d6d7d8", "#f6dbf1", "#a1c0c7",
                                   "#e8ddaf", "#d2d2d2", "#718991", "#b8b3d6", "#a9cbdf",
                                   "#f5dfe1", "#aec6da")))+
  scale_y_continuous(expand = c(0, 0)) +
  xlab("")+
  ylab("")+
  ggtitle("HOS")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        plot.background = element_blank(), 
        panel.border = element_blank(), 
        legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), 
        legend.position = "none")

p_res <- ggplot(result_df_res) +
  geom_area(aes(x, combined_col_res, fill = group),
            position = "fill") +
  scale_fill_manual(values = rev(c("#f4c0bd", "#ead1d1", "#d1eee9", "#f7f4d7", "#efd094",
                                   "#bfdfd2", "#e3edae", "#ccdef4", "#f6e2eb", "#7fd8cc",
                                   "#f4ecb4", "#c2e9e2", "#d6d7d8", "#f6dbf1", "#a1c0c7",
                                   "#e8ddaf", "#d2d2d2", "#718991", "#b8b3d6", "#a9cbdf",
                                   "#f5dfe1", "#aec6da")))+
  scale_y_continuous(expand = c(0, 0)) +
  xlab("")+
  ylab("")+
  ggtitle("RES")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        plot.background = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), 
        legend.position = "none")

p_wtp <- ggplot(result_df_wtp) +
  geom_area(aes(x, combined_col_wtp, fill = group),
            position = "fill") +
  scale_fill_manual(values = rev(c("#f4c0bd", "#ead1d1", "#d1eee9", "#f7f4d7", "#efd094",
                                   "#bfdfd2", "#e3edae", "#ccdef4", "#f6e2eb", "#7fd8cc",
                                   "#f4ecb4", "#c2e9e2", "#d6d7d8", "#f6dbf1", "#a1c0c7",
                                   "#e8ddaf", "#d2d2d2", "#718991", "#b8b3d6", "#a9cbdf",
                                   "#f5dfe1", "#aec6da")))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("")+
  ylab("")+
  ggtitle("WTP")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        plot.background = element_blank(), 
        panel.border = element_blank(), 
        legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), 
        legend.position = "none")

combined_plot <- (p_fli + p_hos + p_res) / p_wtp
print(combined_plot)

ggsave("streamgraph.pdf",combined_plot,device = "pdf",width = 12,height = 6.5)