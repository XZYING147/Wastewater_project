rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure1/ESKAPE")
library(pheatmap)
library(ggplot2)

matrix<-read.table("merged_rpkm_ESKAPE_input.csv",header = T,sep = ",",row.names = 1) 
matrix[is.na(matrix)] <- 0
matrix_log <- log10(matrix+1)

matrix_all<-read.table("merged_rpkm_ESKAPE_input_all.csv",header = T,sep = ",",row.names = 1) 
matrix_all[is.na(matrix_all)] <- 0
matrix_all_log <- log10(matrix_all+1)

eskape_heatmap <- pheatmap(matrix,
                        color = colorRampPalette(c("white","#b23c49" ))(1000), # color?????Զ?????ɫ"#ebeae6
                        # annotation_col = annotation_cols,
                        # annotation_colors = ann_colors,
                        fontsize_col = 1, 
                        cluster_rows = F,
                        cluster_cols = F,
                        show_rownames =T, 
                        show_colnames = F,
                        fontsize = 10,
                        cellwidth=20,
                        cellheight=15,
                        main = "ESKAPE"
) 

eskape_all_heatmap <- pheatmap(matrix_all_log,
                        color = colorRampPalette(c("white","#b23c49" ))(1000), # color?????Զ?????ɫ"#ebeae6
                        # annotation_col = annotation_cols,
                        # annotation_colors = ann_colors,
                        fontsize_col = 1, 
                        cluster_rows = F,
                        cluster_cols = F,
                        show_rownames =T, 
                        show_colnames = F,
                        fontsize = 10,
                        cellwidth=3.5,
                        cellheight=10,
                        main = "ESKAPE"
) 
ggsave("heatmap_ESKAPE.pdf",eskape_heatmap,device = "pdf",width = 7.5,height = 5)
ggsave("heatmap_ESKAPE_all.pdf",eskape_all_heatmap,device = "pdf",width = 10,height = 2.5)
