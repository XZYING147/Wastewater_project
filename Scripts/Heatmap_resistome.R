rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/FigureS3")

library(readr)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

matrix<-read.table("merged_rpkm_type_input.csv",header = T,sep = ",",row.names = 1) 
matrix[is.na(matrix)] <- 0
matrix_log <- log10(matrix+1)

info<-read_excel("metadata_heatmap.xlsx") 
annotation_cols<- data.frame(Group = info$Group, Month = info$Month, Batch = info$Batch, Area = info$Area)
rownames(annotation_cols) <- colnames(matrix)

Group_color <- c("#ecd8aa","#71c7ea","#f0a8c2","#aed2d2","#9abbd2","#dfcebf","#a87471","white") 
names(Group_color) <- c( "Flight","HOS", "RES", "WTP","Ocean","Swine_farm","Sewage","NA") 

Batch_color <- c("#9998ff","#f4a99b","#bd9aad","white")
names(Batch_color) <- c("B1","B2","B3","NA") 

Month_color <- c("#94baf5","#d6c4e0","#e2e2b6","white") 
names(Month_color) <- c("Feb","Mar","Apr","NA") 

Area_color <- c("#f4c0bd","#ead1d1","#d1eee9","#f7f4d7","#efd094","#aec6da","#e3edae","#e9e9e9","#f6e2eb","#7ddacb","#f4ecb4","#c2e9e2", "#ffc0cb","#87ceeb","white")
names(Area_color) <- c("HB","HL","JM","TA","XA","HC","SM","Xinjiang","Hangzhou","Zhejiang","Shanghai","Nanjing","Zhanjiang","Xiamen","NA") 

ann_colors <- list(Group=Group_color,Month=Month_color,Batch=Batch_color,Area=Area_color)

matrix_2<-data.frame(scale(matrix,center = T))
arg_heatmap <- pheatmap(matrix_log,
                        color = colorRampPalette(c("#ebeae6","#b23c49" ))(30),
                        annotation_col = annotation_cols,
                        annotation_colors = ann_colors,
                        fontsize_col = 1, 
                        cluster_rows = F,
                        cluster_cols = F,
                        show_rownames =T, 
                        show_colnames = F,
                        fontsize = 5,
                        cellwidth=2.5,
                        cellheight=10,
) 
ggsave("heatmap_resitome.pdf",arg_heatmap,device = "pdf",width = 25,height = 15)