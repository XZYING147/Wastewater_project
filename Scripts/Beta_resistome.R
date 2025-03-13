rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure2/Beta")

library(ggplot2)
library(dplyr)
library(vegan)
library(cowplot)
library(ggh4x)

samples_metadata <- read.csv("Sample_metadata.csv", fileEncoding = "GBK") %>%
  mutate(
    Sampling_time = as.Date(Sampling_time, format = "%Y/%m/%d"),
    Sampling_year_month = format(Sampling_time, "%Y-%m")
  )

profile_arg <- read.delim("rpkm.subtype_sample.txt")

transpose_and_set_rowname <- function(df) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
  df <- as.data.frame(t(df))
  return(df)
}

profile_arg_t <- transpose_and_set_rowname(profile_arg)
profile_arg_t_relative <- profile_arg_t / rowSums(profile_arg_t)

col_sums_result <- colSums(profile_arg_t_relative)
profile_arg_t_relative <- profile_arg_t_relative[, col_sums_result > 0]
profile_arg_t_relative <- profile_arg_t_relative[, colSums(profile_arg_t_relative) > 0]

Genera_BrayCurtis <- vegdist(profile_arg_t_relative, method = "bray")
pcoa <- cmdscale(Genera_BrayCurtis, k = 7, eig = T, add =TRUE)
pcoa_eig <- pcoa$eig
barplot(pcoa_eig)
pcoa$point
boxplot(pcoa$point)
sample_site <- data.frame(pcoa$point)[1:2]
head(sample_site)
sample_site$Samples <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1','PCoA2')
head(sample_site)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
pcoa_eig
pcoa1 <- paste('PCoA1(', round(100*pcoa_eig[1],2),'%)')
pcoa2 <- paste('PCoA2(', round(100*pcoa_eig[2],2),'%)')
samples_metadata

beta <- merge(sample_site, samples_metadata, by.x = "Samples", by.y = "Sample_ID")

head(beta)

mycol <- c("#ecd8aa", "#71c7ea","#f0a8c2", "#aed2d2")

p1 <- ggplot(beta, aes(PCoA1, PCoA2, gruop = Sampling_type,color = Sampling_type))+
  stat_centroid(aes(xend = PCoA1, yend = PCoA2, colour = Sampling_type),
                geom = "segment", crop_other = F,
                alpha=0.3,size = 1,show.legend = T) +stat_ellipse(aes(color = Sampling_type), size = 0.75, alpha = 1) +
  geom_point(aes(fill=Sampling_type),size = 3, shape = 21, stroke = 1,show.legend = T,color ="black")+
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  theme_test() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
  ) + 
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) 
p1

boxplot_pcoa1 <- ggplot(beta, aes(y = Sampling_type, x = PCoA1, fill = Sampling_type)) +
  geom_boxplot(show.legend = T) +
  scale_fill_manual(values = mycol) +
  labs(y = NULL, x = pcoa1)+  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position="none")+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(axis.title.x = element_text(colour = "black"), axis.text.x = element_text(colour = "black"))
boxplot_pcoa2 <- ggplot(beta, aes(x = Sampling_type, y = PCoA2, fill = Sampling_type)) +
  geom_boxplot(show.legend = T) +
  scale_fill_manual(values = mycol) +
  labs(x = NULL, y = pcoa2) +  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position="none")+
  theme(axis.title.y = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks = element_blank())  # 隐藏 y 轴标签和刻度

final_plot <- plot_grid(
  boxplot_pcoa2, p1,
  NULL, boxplot_pcoa1,
  ncol = 2, nrow = 2,
  rel_widths = c(0.5,1),
  rel_heights = c(1,0.5)
)

final_plot

ggsave("beta_args.pdf",final_plot,device = "pdf",width = 6.5,height = 6.5)
