rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2")

library(ggplot2)
library(reshape2)
library(vegan)
library(Maaslin2)
library(dplyr)

`%notin%` <- Negate(`%in%`)

samples_metadata <- read.csv("Sample_metadata.csv", fileEncoding = "GBK") %>%
  mutate(
    Sampling_time = as.Date(Sampling_time, format = "%Y/%m/%d"),
    Sampling_year_month = format(Sampling_time, "%Y-%m")
  )

pathways <- read.delim("humann3_pathabundance.tsv")
colnames(pathways)[2:ncol(pathways)] <- gsub("_cleaned_concat", "", colnames(pathways)[2:ncol(pathways)])

pathways <- as.data.frame(t(pathways))
colnames(pathways) <- as.character(pathways[1, ]) 
pathways_filtered <- pathways[-1, ]
pathways_filtered$Sample_ID <- rownames(pathways_filtered)

pathways_filtered <- pathways_filtered[, c(ncol(pathways_filtered), 1:(ncol(pathways_filtered)-1))]

rownames(pathways_filtered) <- NULL

pathways_filtered <- pathways[, colnames(pathways) %notin% c("UNMAPPED", "UNINTEGRATED")]
rownames(samples_metadata) <- samples_metadata$Sample_ID
pathways_filtered[pathways_filtered < 0] <- 0
str(pathways_filtered)
str(samples_metadata) 

pathways_filtered[] <- lapply(pathways_filtered, function(x) as.numeric(as.character(x)))
Maaslin2(
  input_data = pathways_filtered,
  input_metadata = samples_metadata,
  fixed_effects = "Sampling_type",
  output = "D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2/microbiome_pathway",
  random_effects = c("Sampling_year_month", "Sampling_area"),
  reference = "Sampling_type,RES",
  normalization = "TSS",
  transform = "LOG"
)

pathways_sign <- read.delim("D:/Desktop/Projects/wastewater/Figure/Figure1/MaAsLin2/microbiome_pathway/significant_results.tsv", header = TRUE, "\t")
pathways_sign$shape <- ifelse(pathways_sign$coef > 0, 1, -1)
pathways_sign <- subset(pathways_sign, qval < 0.05)

p_pathway <- ggplot(pathways_sign, aes(x=coef, y=reorder(feature, coef),color=value)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) + 
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 4)) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway") +
  scale_color_manual(values = c("Flight" = "#ecd8aa", "HOS" ="#71c7ea", "WTP" = "#aed2d2"))

ggsave("pathway_overallsign.pdf", p_pathway, bg = "transparent", width = 15, height = 15)

pathways_sign_WTP <- subset(pathways_sign, value == "WTP")
p_WTP <- ggplot(pathways_sign_WTP, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#aed2d2", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#aed2d2", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_WTP
pathways_sign_HOS <- subset(pathways_sign, value == "HOS")
p_HOS <- ggplot(pathways_sign_HOS, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#71c7ea", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#71c7ea", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_HOS
pathways_sign_Flight <- subset(pathways_sign, value == "Flight")
p_flight <- ggplot(pathways_sign_Flight, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#ecd8aa", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin=coef - stderr, height=0), color = "#ecd8aa", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size=20, face="bold"), axis.text.y = element_text(size=14)) +
  geom_vline(xintercept=0, colour='black', size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p_flight

ggsave("pathway_p_WTP.pdf", p_WTP, bg = "transparent", width = 15, height = 15)
ggsave("pathway_p_HOS.pdf", p_HOS, bg = "transparent", width = 15, height = 15)
ggsave("pathway_p_flight.pdf", p_flight, bg = "transparent", width = 15, height = 15)
