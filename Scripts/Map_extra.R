rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure1/Map")

library(ggplot2)

type_df <- read.csv("location_type.csv")
p1 <- ggplot(data = type_df, aes(Area,value,fill=Sampling_type))+
  geom_bar(stat="identity",position="stack",width=0.8, color = "black", size = 0.05)+
  theme_linedraw() + 
  ylab("Number of samping site") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15),  # x轴标签45度倾斜并末端对齐
    legend.position = "right", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_x_discrete(labels = c("Flight","Haicang", "Jimei", "Tong'an", "Huli", "Siming", "Xiang'an")) +
  scale_color_manual(values = c("#ecd8aa", "#71c7ea", "#f0a8c2", "#aed2d2")) +
  scale_fill_manual(values = c("#ecd8aa", "#71c7ea", "#f0a8c2", "#aed2d2")) +
  xlab(element_blank()) 
p1
ggsave("Map_extra_1.pdf",p1,device = "pdf",width = 10,height = 7.5)

time_df <- read.csv("location_time.csv")
p2 <- ggplot(data = time_df, aes(Area,value,fill=Sampling_time))+
  geom_bar(stat="identity",position="stack",width=0.8, color = "black", size = 0.05)+
  theme_linedraw() + 
  ylab("Number of samping time") +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15),  # x轴标签45度倾斜并末端对齐
    legend.position = "right", 
    axis.title.y = element_text(face = "bold", size = 20)
  ) +
  scale_x_discrete(labels = c("Flight", "Haicang", "Huli", "Jimei", "Siming",  "Tong'an", "Xiang'an")) +
  scale_color_manual(values = c("#e2e2b6","#94baf5","#d6c4e0")) +
  scale_fill_manual(values = c("#e2e2b6","#94baf5","#d6c4e0")) +
  xlab(element_blank()) 
p2
ggsave("Map_extra_2.pdf",p2,device = "pdf",width = 10,height = 7.5)