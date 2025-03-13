rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure3/Chord")

library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(dplyr)
library(tidyr)

data_1 <- read.csv("111.csv")
group_file <- read.csv("type_group.csv")
group <- structure(group_file$group, names = group_file$name)

grid.col <- structure(
  c("WTP_f__Arcobacteraceae"="#eeede7", "RES_f__Arcobacteraceae"="#eeede7", "WTP_f__Burkholderiaceae"="#a4c0bc","HOS_f__Burkholderiaceae_B"="#92b9d4", 
    "WTP_f__Burkholderiaceae_B"="#92b9d4", "WTP_f__Burkholderiaceae_C"="#d9c2de", "FLI_f__Burkholderiaceae_C"="#d9c2de", "WTP_f__Competibacteraceae"="#bd7895", 
    "RES_f__Flavobacteriaceae"="#f6e6eb", "RES_f__Pseudomonadaceae"="#f29f9a", "FLI_f__Pseudomonadaceae"="#f29f9a", "WTP_f__Rhodocyclaceae"="#f9bf83", 
    "HOS_f__Rhodocyclaceae"="#f9bf83", "FLI_f__Ruminococcaceae"="#fbdfc1", "HOS_f__Sphingomonadaceae"="#eddb8f", "RES_f__Sulfurospirillaceae"="#e6b598",
    "HOS_Others"="#cecfd1", "WTP_Others"="#cecfd1", "FLI_Others"="#cecfd1"
    ))

chordDiagram(
  data_1, group = group, grid.col = grid.col,
  col = data_1$Color,
  annotationTrack = c("grid"),
  preAllocateTracks = list(
    track.height = mm_h(4),
    track.margin = c(mm_h(1), 0)
  )
)

circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  for (i in seq_along(group)) {
    if (group[i] == sector.index) {
      circos.points(CELL_META$xcenter, CELL_META$ylim[1] + 0.1, pch = 19, cex = 0.8, col = data_1$Color[i])
    }
  }
})

circos.clear()


data_3 <- read.csv("333.csv")
group_file <- read.csv("type_group.csv")
group <- structure(group_file$group, names = group_file$name)

grid.col_3 <- structure(
  c("WTP_f__Arcobacteraceae"="#eeede7", "RES_f__Arcobacteraceae"="#eeede7", "WTP_f__Burkholderiaceae"="#a4c0bc","HOS_f__Burkholderiaceae_B"="#92b9d4", 
    "WTP_f__Burkholderiaceae_B"="#92b9d4", "WTP_f__Burkholderiaceae_C"="#d9c2de", "FLI_f__Burkholderiaceae_C"="#d9c2de", "WTP_f__Competibacteraceae"="#bd7895", 
    "RES_f__Flavobacteriaceae"="#f6e6eb", "RES_f__Pseudomonadaceae"="#f29f9a", "FLI_f__Pseudomonadaceae"="#f29f9a", "WTP_f__Rhodocyclaceae"="#f9bf83", 
    "HOS_f__Rhodocyclaceae"="#f9bf83", "FLI_f__Ruminococcaceae"="#fbdfc1", "HOS_f__Sphingomonadaceae"="#eddb8f", "RES_f__Sulfurospirillaceae"="#e6b598",
    "HOS_Others"="#cecfd1", "WTP_Others"="#cecfd1", "FLI_Others"="#cecfd1"
  ))

chordDiagram(
  data_3, group = group, grid.col = grid.col_3,
  col = data_3$Color,
  annotationTrack = c("grid"),
  preAllocateTracks = list(
    track.height = mm_h(4),
    track.margin = c(mm_h(1), 0)
  )
)
