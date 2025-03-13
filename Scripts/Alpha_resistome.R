rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure2/Alpha")

library(ggplot2)
library(dplyr)
library(vegan)
library(stringr)
library(ggpubr)
library(cowplot)

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
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)

          },

          draw_group = function(data, panel_scales, coord) {
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))

            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))

            newdata <- rbind(newdata, newdata[1,])

            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },

          draw_key = draw_key_polygon,

          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),

          required_aes = c("x", "y")
  )

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

richness <- as.data.frame(specnumber(profile_arg_t_relative))
richness$Sample <- row.names(richness)
colnames(richness)[1] <-"Richness"
write.csv(richness, "resistome_Richness.csv", row.names = FALSE, quote = TRUE)

shannon_div <- as.data.frame(diversity(profile_arg_t_relative, index = "shannon"))
shannon_div$Sample <- row.names(shannon_div)
colnames(shannon_div)[1] <- "ShannonDiv"
write.csv(shannon_div, "resistome_ShannonDiv.csv", row.names = FALSE, quote = TRUE)

alpha_div <- merge(richness, shannon_div, by = "Sample")
alpha_div <- merge(alpha_div, samples_metadata, by.x = "Sample", by.y = "Sample_ID")
alpha_div_overall <- alpha_div %>%
  mutate(Sampling_year_month = "Overall")
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
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8, size = 1) +
  facet_wrap(~ Sampling_year_month, nrow = 1, strip.position = "top") + 
  theme_linedraw() + 
  ylab("Richness (genera)") +
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
  scale_color_manual(values = c("#f0a8c2", "#71c7ea", "#aed2d2", "#ecd8aa")) + 
  scale_fill_manual(values = c("#f0a8c2", "#71c7ea", "#aed2d2", "#ecd8aa")) + 
  xlab(element_blank()) +
  scale_y_continuous(limits = c(0,900)) + 
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
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8, size = 1) +
  facet_wrap(~ Sampling_year_month, nrow = 1, strip.position = "top") + 
  theme_linedraw() + 
  ylab("shannon(genera)") +
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
  scale_color_manual(values = c("#f0a8c2", "#71c7ea", "#aed2d2", "#ecd8aa")) + 
  scale_fill_manual(values = c("#f0a8c2", "#71c7ea", "#aed2d2", "#ecd8aa")) + 
  xlab(element_blank()) +
  scale_y_continuous(limits = c(0,8)) + 
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
alpha_final
ggsave("alpha_args.pdf",alpha_final,device = "pdf",width = 25,height = 15)
