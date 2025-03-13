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

profile_bact_f <- subset(profile_bact, grepl("^f_", profile_bact$min_level))[, -ncol(profile_bact)]

transpose_and_set_rowname <- function(df) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
  df <- as.data.frame(t(df))
  return(df)
}

profile_bact_f_t <- transpose_and_set_rowname(profile_bact_f)

profile_bact_f_t_relative <- profile_bact_f_t / rowSums(profile_bact_f_t)

col_sums_result <- colSums(profile_bact_f_t_relative)
profile_bact_f_t_relative <- profile_bact_f_t_relative[, col_sums_result > 0]
profile_bact_f_t_relative <- profile_bact_f_t_relative[, colSums(profile_bact_f_t_relative) > 0]

otu.pairwise.adonis <- pairwise.adonis(x=profile_bact_f_t_relative, 
                                       factors=samples_metadata$Sampling_type,
                                       sim.function = "vegdist",
                                       sim.method = "bray",
                                       p.adjust.m = "BH",
                                       reduce = NULL,
                                       perm = 999)

otu.pairwise.adonis <- data.frame(otu.pairwise.adonis, stringsAsFactors = FALSE)

write.csv(otu.pairwise.adonis, 'family_all.pairwise.adonis.csv', row.names = FALSE, quote = FALSE, )
rownames(samples_metadata) <- samples_metadata$Sample_ID

Maaslin2(
  input_data = profile_bact_f_t_relative,
  input_metadata = samples_metadata,
  fixed_effects = "Sampling_type",
  output = "D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2/microbiome_family",
  random_effects = c("Sampling_year_month", "Sampling_area"),
  reference = "Sampling_type,RES",
  normalization = "TSS",
  transform = "LOG"
)
family_sign <- read.delim("D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2/microbiome_family/significant_results.tsv", header = TRUE, "\t")
family_sign$shape <- ifelse(family_sign$coef > 0, 1, -1)
family_sign <- subset(family_sign, qval < 0.05)
p_family <- ggplot(family_sign, aes(x=coef, y=reorder(feature, coef),color=value)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) + #color = "#18678d", 
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 4)) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Genera") +
  scale_color_manual(values = c("Flight" = "#ecd8aa", "HOS" ="#71c7ea", "WTP" = "#aed2d2"))

ggsave("family_overallsign_nm.pdf", p_family, bg = "transparent", width = 10, height = 15)

family_sign_WTP <- subset(family_sign, value == "WTP")
p_WTP <- ggplot(family_sign_WTP, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#aed2d2", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#aed2d2", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_WTP
family_sign_HOS <- subset(family_sign, value == "HOS")
p_HOS <- ggplot(family_sign_HOS, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#71c7ea", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#71c7ea", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_HOS
family_sign_Flight <- subset(family_sign, value == "Flight")
p_flight <- ggplot(family_sign_Flight, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#ecd8aa", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#ecd8aa", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_flight

ggsave("family_p_WTP.pdf", p_WTP, bg = "transparent", width = 15, height = 15)
ggsave("family_p_HOS.pdf", p_HOS, bg = "transparent", width = 15, height = 15)
ggsave("family_p_flight.pdf", p_flight, bg = "transparent", width = 15, height = 15)