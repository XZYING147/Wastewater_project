rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure5/compare")

library(ggplot2)
library(ggpubr)
library(stats)
library(dplyr)

df_res=read.csv("WTP_HOS_vs_RES.txt",  skip = 2, header = TRUE)
p_res <- ggplot() + geom_density(
  data = df_res, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 372, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the WTP-RES vs. HOS-RES lineage-sharing events independent of site", subtitle = "z-score = -4.75 P = 2.62e-6")
ggsave("pic/res.pdf",p_res,device = "pdf",width = 5,height = 5)

df_hos=read.csv("WTP_RES_vs_HOS.txt",  skip = 2, header = TRUE)
p_hos <- ggplot() + geom_density(
  data = df_hos, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 48, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the WTP-HOS vs. RES-HOS lineage-sharing events independent of site", subtitle = "z-score = -6.23 P = 4.61e-10")
ggsave("pic/hos.pdf",p_hos,device = "pdf",width = 5,height = 5)

df_wtp=read.csv("RES_HOS_vs_WTP.txt",  skip = 2, header = TRUE)
p_wtp <- ggplot() + geom_density(
  data = df_wtp, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 372, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the RES-WTP vs. HOS-WTP lineage-sharing events independent of site", subtitle = "z-score = 3.23 P = 0.0012")
ggsave("pic/wtp.pdf",p_wtp,device = "pdf",width = 5,height = 5)

df_all <- read.csv("Different_same.txt",  skip = 2, header = TRUE)
p_all <- ggplot() + geom_density(
  data = df_all, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 802, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the Different site vs. Same site lineage-sharing events", subtitle = "z-score = -5.1 P = 2.66e-7")
ggsave("pic/diff_same.pdf",p_all,device = "pdf",width = 5,height = 5)

df_res_diff <- read.csv("WTP_HOS_vs_RES_diffsite.txt",  skip = 2, header = TRUE)
p_res_diff <- ggplot() + geom_density(
  data = df_res_diff, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 278, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the WTP-RES vs. HOS-RES lineage-sharing events in different site", subtitle = "z-score = -5.9 P = 3.09e-9")
ggsave("pic/res_diff.pdf",p_res_diff,device = "pdf",width = 5,height = 5)

df_hos_diff <- read.csv("WTP_RES_vs_HOS_diffsite.txt",  skip = 2, header = TRUE)
p_hos_diff <- ggplot() + geom_density(
  data = df_hos_diff, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 42, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the WTP-HOS vs. RES-HOS lineage-sharing events in different site", subtitle = "z-score = -6.36 P = 2.9e-10")
ggsave("pic/hos_diff.pdf",p_hos_diff,device = "pdf",width = 5,height = 5)

df_wtp_diff <- read.csv("RES_HOS_vs_WTP_diffsite.txt",  skip = 2, header = TRUE)
p_wtp_diff <- ggplot() + geom_density(
  data = df_wtp_diff, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 278, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the RES-WTP vs. HOS-WTP lineage-sharing events in different site", subtitle = "z-score = 2.05 P = 0.04")
ggsave("pic/wtp_diff.pdf",p_wtp_diff,device = "pdf",width = 5,height = 5)

df_res_same <- read.csv("WTP_HOS_vs_RES_samesite.txt",  skip = 2, header = TRUE)
p_res_same <- ggplot() + geom_density(
  data = df_res_same, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 94, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the WTP-RES vs. HOS-RES lineage-sharing events in same site", subtitle = "z-score = 1.17 P = 0.25")
ggsave("pic/res_same.pdf",p_res_same,device = "pdf",width = 5,height = 5)

df_hos_same <- read.csv("WTP_RES_vs_HOS_samesite.txt",  skip = 2, header = TRUE)
p_hos_same <- ggplot() + geom_density(
  data = df_hos_same, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 6, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the WTP-HOS vs. RES-HOS lineage-sharing events in same site", subtitle = "z-score = -0.45 P = 0.65")
ggsave("pic/hos_same.pdf",p_hos_same,device = "pdf",width = 5,height = 5)

df_wtp_same <- read.csv("RES_HOS_vs_WTP_samesite.txt",  skip = 2, header = TRUE)
p_wtp_same <- ggplot() + geom_density(
  data = df_wtp_same, 
  aes(x = Count), 
  fill = "#18678d",
  alpha = 0.5,
  bw = 1
) + 
  geom_vline(xintercept = 94, color = "black", size = 1.5) + theme_linedraw() + xlab("Lineage-sharing events") +
  labs(title = "Permutation (10,000) of the RES-WTP vs. HOS-WTP lineage-sharing events in same site", subtitle = "z-score = 2.89 P = 0.0034")
ggsave("pic/wtp_same.pdf",p_wtp_same,device = "pdf",width = 5,height = 5)
