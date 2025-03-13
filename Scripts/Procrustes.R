rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/FigureS7")

library(ggplot2)
library(dplyr)
library(vegan)
library(cowplot)
library(ggh4x)
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
Genera_BrayCurtis <- vegdist(profile_bact_g_t_relative, method = "bray")
pcoa <- cmdscale(Genera_BrayCurtis, k = 2, eig = T, add =TRUE)

profile_arg <- read.delim("rpkm.subtype_sample.txt")
profile_arg_t <- transpose_and_set_rowname(profile_arg)
profile_arg_t_relative <- profile_arg_t / rowSums(profile_arg_t)
col_sums_result_arg <- colSums(profile_arg_t_relative)
profile_arg_t_relative <- profile_arg_t_relative[, col_sums_result_arg > 0]
profile_arg_t_relative <- profile_arg_t_relative[, colSums(profile_arg_t_relative) > 0]
Genera_BrayCurtis_arg <- vegdist(profile_arg_t_relative, method = "bray")
pcoa_arg <- cmdscale(Genera_BrayCurtis_arg, k = 2, eig = T, add =TRUE)


pro.s.e <- procrustes(pcoa,pcoa_arg, symmetric = TRUE)
summary(pro.s.e)
plot(pro.s.e, kind = 2)
set.seed(1)
pro.s.e_t <- protest(pcoa,pcoa_arg, permutations = 999)
pro.s.e_t$ss
pro.s.e_t$signif

Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
Pro_X <- data.frame(pro.s.e$rotation)

Pro_Y$Sample <- rownames(Pro_Y)
Pro_Y_metadata <- merge(Pro_Y, samples_metadata, by.x = "Sample", by.y = "Sample_ID")

p_procrustes <- ggplot(data=Pro_Y_metadata) +
  geom_segment(aes(x = X1, y = X2,
                   xend = (X1 + Dim1)/2, yend = (X2 + Dim2)/2),
               size = 0.25) +
  geom_segment(aes(x = (X1 + Dim1)/2, y = (X2 + Dim2)/2,
                   xend = Dim1, yend = Dim2),
               size = 0.25) +
  geom_point(aes(X1, X2, color = Sampling_type), size = 3, shape = 15) +
  geom_point(aes(Dim1, Dim2, color = Sampling_type), size = 3, shape = 16) +
  scale_color_manual(values = c("#f0a8c2", "#71c7ea", "#aed2d2", "#ecd8aa"), limits = c("RES", "HOS", "WTP", "Flight")) +
  theme(panel.grid = element_blank(), # 绘制背景
        panel.background = element_rect(color = 'black',
                                        fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between microbiome and ARG") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[1,2]/Pro_X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[2,2]/Pro_X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\nM2 = 0.7807\nP < 0.001',
           x = -0.075, y = 0.1, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=14,colour = "black",
                                  hjust = 0.5,face = "bold"))
p_procrustes

ggsave("procrustes.pdf",p_procrustes,device = "pdf",width = 8,height = 8)
