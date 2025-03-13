rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/FigureS8")

library(igraph)
library(Hmisc)
library(psych)
library(dplyr)
library(tidyr)

fli <- read.csv("flight_arg_subtype.csv", sep=",", header=T, check.names=F,row.names = 1)
wtp <- read.csv("wtp_arg_subtype.csv", sep=",", header=T, check.names=F,row.names = 1)
res <- read.csv("res_arg_subtype.csv", sep=",", header=T, check.names=F,row.names = 1)
hos <- read.csv("hos_arg_subtype.csv", sep=",", header=T, check.names=F,row.names = 1)
group <- read.table("group_outer.txt", sep="\t", header=T, check.names=F)

fli <- as.data.frame(t(fli))
fli$sample <- rownames(fli)
wtp <- as.data.frame(t(wtp))
wtp$sample <- rownames(wtp)
res <- as.data.frame(t(res))
res$sample <- rownames(res)
hos <- as.data.frame(t(hos))
hos$sample <- rownames(hos)

df_fli_hos <- merge(fli, hos, by = "sample")
df_fli_hos_res <- merge(df_fli_hos, res, by = "sample")
df <- merge(df_fli_hos_res, wtp, by = "sample")
rownames(df) <- df$sample
df <- df[-1]
head(df)

data<-as.matrix(df)
cor<- corr.test(data, method="spearman",adjust="BH")
data.cor <- as.data.frame(cor$r)

r.cor<-data.frame(cor$r)
p.cor<-data.frame(cor$p)
r.cor[p.cor>0.05] <- 0
r.cor[abs(r.cor) < 0.7] <- 0
r.cor[r.cor == 1] <- 0
r.cor[r.cor == "1"] <- 0

r.cor$from = rownames(r.cor)
p.cor$from = rownames(p.cor)
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
write.csv(cor.data, "cor_outer_0.7.csv")

vertices <- c(as.character(cor.data$from),as.character(cor.data$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarise()
colnames(vertices) <- "name"
vertices <- vertices %>%
  left_join(group,by="name")
vertices$group <- factor(vertices$group, levels = c("Flight","WTP","RES","HOS"))
vertices <- vertices %>%
  arrange(group)

graph <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
E(graph)$weight <- abs(E(graph)$r)
V(graph)$label <- V(graph)$name
write_graph(graph, "cor_outer_0.7.graphml", format="graphml")
