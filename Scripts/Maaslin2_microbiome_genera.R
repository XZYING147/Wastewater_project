rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2")

library(ggplot2)
library(stringr)
library(dplyr)
library(vegan)
library(Maaslin2)
library(pairwiseAdonis)

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

otu.pairwise.adonis <- pairwise.adonis(x=profile_bact_g_t_relative, 
                                       factors=samples_metadata$Sampling_type,
                                       sim.function = "vegdist",
                                       sim.method = "bray",
                                       p.adjust.m = "BH",
                                       reduce = NULL,
                                       perm = 999)

otu.pairwise.adonis <- data.frame(otu.pairwise.adonis, stringsAsFactors = FALSE)

write.csv(otu.pairwise.adonis, 'genera_sall.pairwise.adonis.csv', row.names = FALSE, quote = FALSE, )
rownames(samples_metadata) <- samples_metadata$Sample_ID

Maaslin2(
  input_data = profile_bact_g_t_relative,
  input_metadata = samples_metadata,
  fixed_effects = "Sampling_type",
  output = "D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2/microbiome_genera",
  random_effects = c("Sampling_year_month", "Sampling_area"),
  reference = "Sampling_type,RES",
  normalization = "TSS",
  transform = "LOG"
)
genera_sign <- read.delim("D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2/microbiome_genera/significant_results.tsv", header = TRUE, "\t")
genera_sign$shape <- ifelse(genera_sign$coef > 0, 1, -1)
genera_sign <- subset(genera_sign, qval < 0.05)
p_test_nm <- ggplot(genera_sign, aes(x=coef, y=reorder(feature, coef),color=value)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) + #color = "#18678d", 
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 4)) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Genera") +
  scale_color_manual(values = c("Flight" = "#ecd8aa", "HOS" ="#71c7ea", "WTP" = "#aed2d2"))

ggsave("genera_overallsign_nm.pdf", p_test_nm, bg = "transparent", width = 10, height = 15)

genera_sign_WTP <- subset(genera_sign, value == "WTP")
p_WTP <- ggplot(genera_sign_WTP, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#aed2d2", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#aed2d2", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_WTP
genera_sign_HOS <- subset(genera_sign, value == "HOS")
p_HOS <- ggplot(genera_sign_HOS, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#71c7ea", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#71c7ea", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_HOS
genera_sign_Flight <- subset(genera_sign, value == "Flight")
p_flight <- ggplot(genera_sign_Flight, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#ecd8aa", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#ecd8aa", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_flight

ggsave("genera_p_WTP.pdf", p_WTP, bg = "transparent", width = 15, height = 15)
ggsave("genera_p_HOS.pdf", p_HOS, bg = "transparent", width = 15, height = 15)
ggsave("genera_p_flight.pdf", p_flight, bg = "transparent", width = 15, height = 15)