rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure2/AVD")

library(ggplot2)
library(ggdist)
library(gghalves)

otu1 <- read.delim("rpkm.subtype_sample_AVD.txt",row.names = 1)

ai <- abs(otu1 - apply(otu1, 1, mean)) / apply(otu1, 1, sd)

ai[is.na(ai)] <- 0

avd <- colSums(ai, na.rm = TRUE) / (1 * nrow(otu1))
group <- read.csv("annotation_4_AVD.csv", check.names = FALSE,stringsAsFactors = T)
group$AVD <- avd
group_ALL<- group[-3,]

avd_all<- ggplot(data = group_ALL, mapping = aes(x = Process, y = AVD, fill = Process)) +
  stat_gradientinterval(position = "dodge",colour = NA, width = 0.8, fill_type = "segments") +
  stat_halfeye(adjust = .3,width = .3,.width = 0,
               justification = -.3,point_colour = 'NA',
               slab_fill=NA, slab_colour='#3e2c12',
               slab_size=0.4)+
  geom_boxplot(width = .15,outlier.shape = NA,fill=c("#ecd8aa", "#71c7ea", "#f0a8c2", "#aed2d2")) +
  scale_fill_manual(values=c("#ecd8aa", "#71c7ea", "#f0a8c2", "#aed2d2"))+
  geom_half_point(side  = "l",alpha = 0.5, size  = 1.8,aes(color = Process)) +
  scale_color_manual(values = c("#ecd8aa", "#71c7ea", "#f0a8c2", "#aed2d2")) +
  guides(fill="none", alpha='none') +
  theme_light()+
  theme(panel.grid = element_line(color = 'NA'),
        panel.grid.minor  = element_line(color ='NA'),
        panel.border = element_rect(fill = 'NA', color = 'black', size = 0.8, linetype = 'solid'),
        axis.text.x =element_text(size = 8, angle = 0,vjust = 0.6,color = 'black',hjust = 0.5,
                                  face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, angle = 0,vjust = 0.6,color = 'black',hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(vjust=0.2,size = 13))+ theme(legend.position="none")

avd_all
ggsave("avd_all.pdf",avd_all,device = "pdf",width = 7.5,height = 10)
