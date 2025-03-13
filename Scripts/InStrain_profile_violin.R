rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure5/diversity_snv")

library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)
library(stringr)

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


trans_tsv_files <- list.files(path = "IS_profile_trans", pattern = "\\.tsv$", full.names = TRUE)
if(length(trans_tsv_files) == 0) stop("TSV missing")
trans_combined_df <- data.frame()
for(file in trans_tsv_files){
  tryCatch({
    df <- read_tsv(file, col_types = cols(.default = col_character())) %>%
      type_convert()
    
    required_cols <- c("coverage", "breadth")
    if(!all(required_cols %in% colnames(df))){
      warning(paste("File", basename(file), "Missing, passed!"))
      next
    }
    filtered <- df %>%
      filter(
        as.numeric(coverage) > 10,
        as.numeric(breadth) > 0.8
        ,as.numeric(nucl_diversity) < 0.015
      ) %>%
      mutate(
        source = str_remove(basename(file), "\\.tsv$")
      )
    trans_combined_df <- bind_rows(trans_combined_df, filtered)
  }, error = function(e){
    message(paste(basename(file), e$message))
  })
}

trans_comparisons <- combn(as.character(unique(trans_combined_df$source)), 2, simplify = FALSE)

p_trans_diversity <- ggplot(trans_combined_df, aes(x = source, y = nucl_diversity, color = source)) +
  geom_boxplot(aes(fill = source), color = "black", width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = source), position = position_nudge(x = 0.24), alpha = 0.7) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8, size = 1) +
  theme_linedraw() + 
  ylab("Diversity") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15), 
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_color_manual(values = c("Intra"="#a3ced6","Inter"="#dfaeb6")) + 
  scale_fill_manual(values = c("Intra"="#a3ced6","Inter"="#dfaeb6")) + 
  xlab(element_blank())+
  stat_compare_means(
    comparisons = trans_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    exact = FALSE,
    hide.ns = TRUE
  )
p_trans_diversity
ggsave("trans_diversity.pdf", p_trans_diversity, device = "pdf", width = 4, height = 10)

p_trans_snv <- ggplot(trans_combined_df, aes(x = source, y = SNV_distance_mean, color = source)) +
  geom_boxplot(aes(fill = source), color = "black", width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = source), position = position_nudge(x = 0.24), alpha = 0.7) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8, size = 1) +
  theme_linedraw() + 
  ylab("SNV_distance_mean") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15), 
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_color_manual(values = c("Intra"="#a3ced6","Inter"="#dfaeb6")) + 
  scale_fill_manual(values = c("Intra"="#a3ced6","Inter"="#dfaeb6")) +  
  xlab(element_blank())+
  stat_compare_means(
    comparisons = trans_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    exact = FALSE,
    hide.ns = TRUE
  )
p_trans_snv
ggsave("trans_snv_distance.pdf", p_trans_snv, device = "pdf", width = 4, height = 10)



type_tsv_files <- list.files(path = "IS_profile_type", pattern = "\\.tsv$", full.names = TRUE)
if(length(type_tsv_files) == 0) stop("TSV missing")
type_combined_df <- data.frame()
for(file in type_tsv_files){
  tryCatch({
    df <- read_tsv(file, col_types = cols(.default = col_character())) %>%
      type_convert()
    
    required_cols <- c("coverage", "breadth")
    if(!all(required_cols %in% colnames(df))){
      warning(paste("File", basename(file), "Missing, passed!"))
      next
    }
    filtered <- df %>%
      filter(
        as.numeric(coverage) > 10,
        as.numeric(breadth) > 0.8
        ,as.numeric(nucl_diversity) < 0.015
      ) %>%
      mutate(
        source = str_remove(basename(file), "\\.tsv$")
      )
    type_combined_df <- bind_rows(type_combined_df, filtered)
  }, error = function(e){
    message(paste(basename(file), e$message))
  })
}

type_comparisons <- combn(as.character(unique(type_combined_df$source)), 2, simplify = FALSE)

p_type_diversity <- ggplot(type_combined_df, aes(x = source, y = nucl_diversity, color = source)) +
  geom_boxplot(aes(fill = source), color = "black", width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = source), position = position_nudge(x = 0.24), alpha = 0.7) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8, size = 1) +
  theme_linedraw() + 
  ylab("Diversity") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15), 
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_color_manual(values = c("FLI"="#ecd9ae", "RES"="#f0acc5", "HOS"="#78c9eb", "WTP"="#b2d4d4")) + 
  scale_fill_manual(values = c("FLI"="#ecd9ae", "RES"="#f0acc5", "HOS"="#78c9eb", "WTP"="#b2d4d4")) + 
  xlab(element_blank())+
  stat_compare_means(
    comparisons = type_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    exact = FALSE,
    hide.ns = TRUE
  )
p_type_diversity
ggsave("type_diversity.pdf", p_type_diversity, device = "pdf", width = 7.5, height = 10)

p_type_snv <- ggplot(type_combined_df, aes(x = source, y = SNV_distance_mean, color = source)) +
  geom_boxplot(aes(fill = source), color = "black", width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = source), position = position_nudge(x = 0.24), alpha = 0.7) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8, size = 1) +
  theme_linedraw() + 
  ylab("SNV_distance_mean") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15), 
    legend.position = "none", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_color_manual(values = c("FLI"="#ecd9ae", "RES"="#f0acc5", "HOS"="#78c9eb", "WTP"="#b2d4d4")) + 
  scale_fill_manual(values = c("FLI"="#ecd9ae", "RES"="#f0acc5", "HOS"="#78c9eb", "WTP"="#b2d4d4")) +  
  xlab(element_blank())+
  stat_compare_means(
    comparisons = type_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    exact = FALSE,
    hide.ns = TRUE
  )
p_type_snv
ggsave("type_snv_distance.pdf", p_type_snv, device = "pdf", width = 7.5, height = 10)
