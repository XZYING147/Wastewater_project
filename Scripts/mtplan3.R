setwd("C:/Users/Jobs/Nutstore/2/zhiying_jobs study/wastewater")
setwd("C:/Users/XJB/Desktop/zhiying_jobs study/wastewater/")
`%notin%` <- Negate(`%in%`)
library(ggplot2)
library(reshape2)
library(vegan)
library(Maaslin2)
library(dplyr)
library(stringr)
library(ggpubr)
library(PupillometryR)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

# remotes::install_github("biobakery/maaslin2")
# 读取样本元数据
samples_metadata <- read.csv("wordfile/Sample_metadata.csv", fileEncoding = "GBK") %>%
  mutate(
    Sampling_time = as.Date(Sampling_time, format = "%Y/%m/%d"),
    Sampling_year_month = format(Sampling_time, "%Y-%m")
  )
# 读取路径丰度数据
profile_bact <- read.delim("data/metaphlan4/metaphlan4_taxonomy.tsv")
# 提取 ID 列中最后一个分类级别
profile_bact$min_level <- str_extract(profile_bact$ID, "(k__|p__|c__|o__|f__|g__|s__|t__)[^|]*$")
# 创建各分类级别的数据框
profile_bact_k <- subset(profile_bact, grepl("^k__", profile_bact$min_level))[, -ncol(profile_bact)]
profile_bact_p <- subset(profile_bact, grepl("^p__", profile_bact$min_level))[, -ncol(profile_bact)]
profile_bact_c <- subset(profile_bact, grepl("^c__", profile_bact$min_level))[, -ncol(profile_bact)]
profile_bact_o <- subset(profile_bact, grepl("^o__", profile_bact$min_level))[, -ncol(profile_bact)]
profile_bact_f <- subset(profile_bact, grepl("^f__", profile_bact$min_level))[, -ncol(profile_bact)]
profile_bact_g <- subset(profile_bact, grepl("^g__", profile_bact$min_level))[, -ncol(profile_bact)]
profile_bact_s <- subset(profile_bact, grepl("^s__", profile_bact$min_level))[, -ncol(profile_bact)]
profile_bact_t <- subset(profile_bact, grepl("^t__", profile_bact$min_level))[, -ncol(profile_bact)]

# 定义一个函数，用于转置数据框并将第一列设置为行名
transpose_and_set_rowname <- function(df) {
  # 将第一列设为行名
  rownames(df) <- df[, 1]
  # 删除第一列
  df <- df[, -1]
  # 转置数据框
  df <- as.data.frame(t(df))
  return(df)
}

# 对每个 profile_bact_x 数据框应用该函数
profile_bact_k_t <- transpose_and_set_rowname(profile_bact_k)
profile_bact_p_t <- transpose_and_set_rowname(profile_bact_p)
profile_bact_c_t <- transpose_and_set_rowname(profile_bact_c)
profile_bact_o_t <- transpose_and_set_rowname(profile_bact_o)
profile_bact_f_t <- transpose_and_set_rowname(profile_bact_f)
profile_bact_g_t <- transpose_and_set_rowname(profile_bact_g)
profile_bact_s_t <- transpose_and_set_rowname(profile_bact_s)
profile_bact_t_t <- transpose_and_set_rowname(profile_bact_t)

profile_bact_g_t_relative <- profile_bact_g_t / rowSums(profile_bact_g_t) %>%
  profile_bact_g_t[, colSums(profile_bact_g_t) > 0]
profile_bact_g_t_relative <- profile_bact_g_t_relative[, colSums(profile_bact_g_t_relative) > 0]


richness <- as.data.frame(specnumber(profile_bact_g_t_relative))
richness$Sample <- row.names(richness)
colnames(richness)[1] <-"Richness"

shannon_div <- as.data.frame(diversity(profile_bact_g_t_relative, index = "shannon"))
shannon_div$Sample <- row.names(shannon_div)
colnames(shannon_div)[1] <- "ShannonDiv"

alpha_div <- merge(richness, shannon_div, by = "Sample")
alpha_div <- merge(alpha_div, samples_metadata, by.x = "Sample", by.y = "Sample_ID")
alpha_div_overall <- alpha_div %>%
  mutate(Sampling_year_month = "Overall")  # 将所有数据的月份标记为“Overall”
alpha_div_extended <- bind_rows(alpha_div, alpha_div_overall)
alpha_div_extended$Richness <- as.numeric(as.character(alpha_div_extended$Richness))
comparisons <- combn(as.character(unique(alpha_div_extended$Sampling_type)), 2, simplify = FALSE)
alpha_div_extended$Sampling_year_month <- factor(alpha_div_extended$Sampling_year_month, 
                                                 levels = c("2023-02", "2023-03", "2023-04", "Overall"),
                                                 labels = c("2023-Feb", "2023-Mar", "2023-Apr", "Overall"))
alpha_div_extended$Sampling_type <- factor(alpha_div_extended$Sampling_type, 
                                                 levels = c( "RES", "HOS", "WTP","Flight"))
p1_alpha_rich <- ggplot(alpha_div_extended, aes(x = Sampling_type, y = Richness, color = Sampling_type)) +
  geom_boxplot(aes(fill = Sampling_type), color = "black", width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Sampling_type), position = position_nudge(x = 0.24), alpha = 0.7) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ Sampling_year_month, nrow = 1, strip.position = "top") + 
  theme_linedraw() + 
  ylab("Richness (genera)") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 5), 
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),  # x轴标签45度倾斜并末端对齐
    strip.background = element_rect(fill = "black"),  # 分面标题背景设为黑色
    strip.text = element_text(color = "white", face = "bold", size = 7.5),  # 分面标题文字设为白色
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 10)
  ) +
  scale_x_discrete(labels = c("Residence", "WWTP", "Hospital", "Flight")) + 
  scale_color_manual(values = c("#873e23", "#18678d", "#0a7c40", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#0a7c40", "#626262")) + 
  xlab(element_blank()) +
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", 
    p.adjust.method = "BH", 
    exact = FALSE
  )
p1_alpha_rich
p1_alpha_shannon<- ggplot(alpha_div_extended, aes(x = Sampling_type, y = ShannonDiv, color = Sampling_type)) +
  geom_boxplot(aes(fill = Sampling_type), color = "black", width = 0.4, size = 1, alpha = 0.5,stroke = 0.5) + 
  geom_flat_violin(aes(fill = Sampling_type), position = position_nudge(x = 0.24), alpha = 0.7) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ Sampling_year_month, nrow = 1, strip.position = "top") + 
  theme_linedraw() + 
  ylab("shannon(genera)") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 5), 
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),  # x轴标签45度倾斜并末端对齐
    strip.background = element_rect(fill = "black"),  # 分面标题背景设为黑色
    strip.text = element_text(color = "white", face = "bold", size = 7.5),  # 分面标题文字设为白色
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 5)
  ) +
  scale_x_discrete(labels = c("Residence", "WWTP", "Hospital", "Flight")) + 
  scale_color_manual(values = c("#873e23", "#18678d", "#0a7c40", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#0a7c40", "#626262")) + 
  xlab(element_blank()) +
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", 
    p.adjust.method = "BH", 
    exact = FALSE
  )
p1_alpha_shannon
alpha_final <- plot_grid(
  p1_alpha_rich, p1_alpha_shannon,
  ncol = 2, nrow = 1)

ggsave("alpha.pdf",alpha_final,device = "pdf",width = 10,height = 4)

Genera_BrayCurtis <- vegdist(profile_bact_g_t_relative, method = "bray")#计算距离，需要确保每个library都是非0值
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

mycol <- c("#873e23", "#18678d", "#0a7c40", "#626262")
library(ggh4x)
p2 <- ggplot(beta, aes(PCoA1, PCoA2, gruop = Sampling_type,color = Sampling_type))+
  stat_centroid(aes(xend = PCoA1, yend = PCoA2, colour = Sampling_type),
                geom = "segment", crop_other = F,
                alpha=0.3,size = 1,show.legend = F) +stat_ellipse(aes(color = Sampling_type), size = 0.75, alpha = 1) +
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
p2

boxplot_pcoa1 <- ggplot(beta, aes(y = Sampling_type, x = PCoA1, fill = Sampling_type)) +
  geom_boxplot(show.legend = F) +
  scale_fill_manual(values = mycol) +
  labs(y = NULL, x = pcoa1)+  theme_classic2() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position="none")+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(axis.title.x = element_text(colour = "black"), axis.text.x = element_text(colour = "black"))
boxplot_pcoa2 <- ggplot(beta, aes(x = Sampling_type, y = PCoA2, fill = Sampling_type)) +
  geom_boxplot(show.legend = F) +
  scale_fill_manual(values = mycol) +
  labs(x = NULL, y = pcoa2) +  theme_classic2() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position="none")+
  theme(axis.title.y = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks = element_blank())  # 隐藏 y 轴标签和刻度


library(cowplot)
# 将图形排列在一起
final_plot <- plot_grid(
  boxplot_pcoa2, p2,
  NULL, boxplot_pcoa1,
  ncol = 2, nrow = 2,
  rel_widths = c(0.5,1),
  rel_heights = c(1,0.5)
)

final_plot

ggsave("beta.pdf",final_plot,device = "pdf",width = 6.5,height = 6.5)

otu.pairwise.adonis <- pairwise.adonis(x=stat_compare_means, factors=samples_metadata$Sampling_type,
                                       sim.function = "vegdist",
                                       sim.method = "bray",
                                       p.adjust.m = "BH",
                                       reduce = NULL,
                                       perm = 999)

otu.pairwise.adonis <- data.frame(otu.pairwise.adonis, stringsAsFactors = FALSE)

write.csv(otu.pairwise.adonis, 'all.pairwise.adonis.csv', row.names = FALSE, quote = FALSE, )
rownames(samples_metadata) <- samples_metadata$Sample_ID


Maaslin2(
  input_data = profile_bact_g_t_relative,
  input_metadata = samples_metadata,
  fixed_effects = "Sampling_type",
  output = "C:/Study1/zhiying",
  random_effects = c("Sampling_year_month", "Sampling_area"),
  reference = "Sampling_type,RES",
  normalization = "TSS",
  transform = "LOG"
)
pathways_sign <- read.delim("C:/Study1/zhiying/significant_results.tsv", header = TRUE, "\t")
pathways_sign$shape <- ifelse(pathways_sign$coef > 0, 1, -1)
pathways_sign <- subset(pathways_sign, qval < 0.05)
p_test_nm <- ggplot(pathways_sign, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_test <- ggplot(pathways_sign, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#18678d", size = 0.7) +
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_blank()) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted")+
  scale_y_discrete("")


ggsave("genera_overallsign_nm.pdf", p_test_nm, bg = "transparent", width = 15, height = 15)
ggsave("genera_overallsign.pdf", p_test, bg = "transparent", width = 15, height = 15)

pathways_sign_WTP <- subset(pathways_sign, value == "WTP")
p_WTP <- ggplot(pathways_sign_WTP, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_WTP
pathways_sign_HOS <- subset(pathways_sign, value == "HOS")
p_HOS <- ggplot(pathways_sign_HOS, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_HOS
pathways_sign_Flight <- subset(pathways_sign, value == "Flight")
p_flight <- ggplot(pathways_sign_Flight, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_flight

ggsave("genera_p_WTP.pdf", p_WTP, bg = "transparent", width = 15, height = 15)
ggsave("genera_p_HOS.pdf", p_HOS, bg = "transparent", width = 15, height = 15)
ggsave("genera_p_flight.pdf", p_flight, bg = "transparent", width = 15, height = 15)


