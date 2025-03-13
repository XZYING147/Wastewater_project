rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure2/Network")

library(igraph)
library(Hmisc)
library(psych)
library(dplyr)
library(tidyr)

network <- read.csv("rpkm.subtype_sample.csv", sep=",", header=T,row.names = 1)
group <- read.csv("group_subtype_sample_part.csv", sep=",", header=T)

data<-as.matrix(network)

cor<- corr.test(data, method="spearman",adjust="BH")

data.cor <- as.data.frame(cor$r)

r.cor<-data.frame(cor$r)
p.cor<-data.frame(cor$p)
r.cor[p.cor>0.05] <- 0
r.cor[abs(r.cor) < 0.7] <- 0
r.cor[r.cor == 1] <- 0
r.cor[r.cor == "1"] <- 0

r.cor$from = rownames(r.cor)[1:1552,1553:1804]
p.cor$from = rownames(p.cor)[1:1552,1553:1804]
p_value <-  p.cor %>%
  gather(key = "to", value = "p", -from) %>%
  data.frame()
p_value <- p_value[, -3]
cor.data<- r.cor %>%
  gather(key = "to", value = "r", -from) %>%
  data.frame() %>%
  left_join(p_value, by=c("from","to")) %>%
  mutate(
    linecolor = ifelse(r > 0,"positive","negative"),
    linesize = abs(r)
  )

cor.data <- cor.data[abs(cor.data$r)>0.7, ]
write.csv(cor.data, "correlations_all_0.7.csv")


filtered_correlations_all_0_5 <- read.csv("correlations_all_0.7.csv", sep=",", header=T,row.names = 1)
cor.data <- filtered_correlations_all_0_5

vertices <- c(as.character(cor.data$from),as.character(cor.data$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarise()
colnames(vertices) <- "name"
vertices <- vertices %>%
  left_join(group,by="name")
vertices$group <- factor(vertices$group, levels = c("Flight","WTP","RES","HOS","aminoglycoside","beta_lactam","chloramphenicol","fosfomycin","MLS","multidrug","puromycin","quinolone","rifamycin","tetracycline","trimethoprim","vancomycin"))
vertices <- vertices %>%
  arrange(group)

graph <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
E(graph)$weight <- abs(E(graph)$r)*3
V(graph)$label <- V(graph)$name
write_graph(graph, "network_test.graphml", format="graphml")
