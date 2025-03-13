rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure4/GLM")

library(jtools)
library(sjPlot)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(piecewiseSEM)
library(MuMIn)
library(glmm.hp)

virus_data <- read.csv("GLM_sewage.csv",check.names = F, header = T, row.names = 1)

glm_model_1 <- glm(Resistome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                     Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                     Latitude + Longitude + Max_Population + Mean_Population + Microbiome_Shannon + pH + Rainfall + 
                     `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                   data = virus_data, family = poisson)
glm_model_2 <- glm(Resistome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                     Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                     Latitude + Longitude + Max_Population + Mean_Population + Microbiome_Shannon + pH + Rainfall + 
                     `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                   data = virus_data, family = Gamma(link = "inverse"))
glm_model_3 <- glm(Resistome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                     Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                     Latitude + Longitude + Max_Population + Mean_Population + Microbiome_Shannon + pH + Rainfall + 
                     `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                   data = virus_data, family = inverse.gaussian(link = "1/mu^2"))
glm_model_4 <- glm(Resistome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                     Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                     Latitude + Longitude + Max_Population + Mean_Population + Microbiome_Shannon + pH + Rainfall + 
                     `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                   data = virus_data, family = quasi(link = "identity", variance = "constant"))
glm_model_5 <- glm(Resistome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                     Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                     Latitude + Longitude + Max_Population + Mean_Population + Microbiome_Shannon + pH + Rainfall + 
                     `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                   data = virus_data, family = quasipoisson(link = "log"))
glm_model_6 <- glm(Resistome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                     Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                     Latitude + Longitude + Max_Population + Mean_Population + Microbiome_Shannon + pH + Rainfall + 
                     `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                   data = virus_data, family = gaussian)

AIC_1 <- AIC(glm_model_1)
AIC_2 <- AIC(glm_model_2)
AIC_3 <- AIC(glm_model_3)
AIC_4 <- AIC(glm_model_4)
AIC_5 <- AIC(glm_model_5)
AIC_6 <- AIC(glm_model_6)

min_AIC <- min(AIC_1, AIC_2, AIC_3, AIC_4, AIC_5, AIC_6)

if (min_AIC == AIC_1) {
  best_model <- "glm_model_1"
} else if (min_AIC == AIC_2) {
  best_model <- "glm_model_2"
} else if (min_AIC == AIC_3) {
  best_model <- "glm_model_3"
} else if (min_AIC == AIC_4) {
  best_model <- "glm_model_4"
} else if (min_AIC == AIC_5) {
  best_model <- "glm_model_5"
} else if (min_AIC == AIC_6) {
  best_model <- "glm_model_6"
}
cat("Best Model:", best_model, "\n")

glm_model_6_conditional_r2 <- rsquared(glm_model_6, method = "nagelkerke")
glm_model_6_marginal_r2 <- r.squaredGLMM(glm_model_6)

p1_all <-plot_model(glm_model_6, type = "slope", colors = c("#c15560","#7fa4d1"),
                    show.intercept = T,show.p = TRUE,  grid.breaks=NULL,
                    show.values = T, value.offset = 5)+theme_light()+
  theme(panel.grid = element_line(color = 'NA'),
        panel.grid.minor  = element_line(color ='NA'),
        panel.border = element_rect(fill = 'NA', color = 'black', size = 0.8, linetype = 'solid'),
        axis.text.x =element_text(size = 16, angle = 0,vjust = 0.6,color = 'black',hjust = 0.5,
                                  face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, angle = 0,vjust = 0.6,color = 'black',hjust = 1, face = 'bold'),
        axis.title.y = element_text(vjust=0.2,size = 13))+   theme(legend.position="none", strip.text = element_text(color = "white", size = 14))
p1_all
ggsave("GLM_resistome_sup.pdf",p1_all,device = "pdf",width = 20,height = 15)


p2_all <-plot_model(glm_model_6, type = "std2", colors = c("#c15560","#7fa4d1"),
                    title = "GLM estimates of effect on resistome",show.intercept = F,show.legend = T,
                    show.p = TRUE,  grid.breaks=NULL,
                    show.values = T, value.offset = 0.35)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_light()+
  theme(panel.grid = element_line(color = 'NA'),
        panel.grid.minor  = element_line(color ='NA'),
        panel.border = element_rect(fill = 'NA', color = 'black', size = 0.8, linetype = 'solid'),
        axis.text.x =element_text(size = 16, angle = 0,vjust = 0.6,color = 'black',hjust = 0.5,
                                  face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, angle = 0,vjust = 0.25,color = 'black',hjust = 1, face = 'bold'),
        axis.title.y = element_text(vjust=0.2,size = 13))+ theme(legend.position="left")

p2_all 
ggsave("GLM_resistome.pdf",p2_all,device = "pdf",width = 12,height = 15)

summary(glm_model_6)
summary_text <- capture.output(summary(glm_model_6))
write.csv(summary_text,file = "GLM_resistome_parameter.csv")



micro_glm_model_1 <- glm(Microbiome_Richness ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                           Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                           Latitude + Longitude + Max_Population + Mean_Population + Resistome_Shannon + pH + Rainfall + 
                           `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                         data = virus_data, family = poisson)

micro_glm_model_2 <- glm(Microbiome_Richness ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                           Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                           Latitude + Longitude + Max_Population + Mean_Population + Resistome_Shannon + pH + Rainfall + 
                           `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                         data = virus_data, family = Gamma(link = "inverse"))

micro_glm_model_3 <- glm(Microbiome_Richness ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                           Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                           Latitude + Longitude + Max_Population + Mean_Population + Resistome_Shannon + pH + Rainfall + 
                           `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                         data = virus_data, family = inverse.gaussian(link = "1/mu^2"))

micro_glm_model_4 <- glm(Microbiome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                           Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                           Latitude + Longitude + Max_Population + Mean_Population + Resistome_Shannon + pH + Rainfall + 
                           `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                         data = virus_data, family = quasi(link = "identity", variance = "constant"))

micro_glm_model_5 <- glm(Microbiome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp + 
                           Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time + 
                           Latitude + Longitude + Max_Population + Mean_Population + Resistome_Shannon + pH + Rainfall + 
                           `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                         data = virus_data, family = quasipoisson(link = "log"))

micro_glm_model_6 <- glm(Microbiome_Shannon ~ Average_Humidity + Average_Water_Temp + Average_Temp +
                           Chemical_Oxygen_Demand + Color_Fold + FlowRate + Hydraulic_Retention_Time +
                           Latitude + Longitude + Max_Population + Mean_Population + Resistome_Shannon + pH + Rainfall +
                           `Total_Ammonia-Nitrogen` + Total_Chlorine + Total_Nitrogen + Total_Phosphorus + Total_Suspended_Solid,
                         data = virus_data, family = gaussian)

micro_AIC_1 <- AIC(micro_glm_model_1)
micro_AIC_2 <- AIC(micro_glm_model_2)
micro_AIC_3 <- AIC(micro_glm_model_3)
micro_AIC_4 <- AIC(micro_glm_model_4)
micro_AIC_5 <- AIC(micro_glm_model_5)
micro_AIC_6 <- AIC(micro_glm_model_6)

micro_min_AIC <- min(micro_AIC_1, micro_AIC_2, micro_AIC_3, micro_AIC_4, micro_AIC_5, micro_AIC_6)

if (micro_min_AIC == micro_AIC_1) {
  micro_best_model <- "glm_model_1"
} else if (min_AIC == micro_AIC_2) {
  micro_best_model <- "glm_model_2"
} else if (min_AIC == micro_AIC_3) {
  micro_best_model <- "glm_model_3"
} else if (min_AIC == micro_AIC_4) {
  micro_best_model <- "glm_model_4"
} else if (min_AIC == micro_AIC_5) {
  micro_best_model <- "glm_model_5"
} else if (min_AIC == micro_AIC_6) {
  micro_best_model <- "glm_model_6"
}
cat("Best Model:", micro_best_model, "\n")

micro_glm_model_6_conditional_r2 <- rsquared(micro_glm_model_6, method = "nagelkerke")
micro_glm_model_6_marginal_r2 <- r.squaredGLMM(micro_glm_model_6)

micro_p1_all <-plot_model(micro_glm_model_6, type = "slope", colors = c("#c15560","#7fa4d1"),
                          show.intercept = T,show.p = TRUE,  grid.breaks=NULL,
                          show.values = T, value.offset = 5)+theme_light()+
  theme(panel.grid = element_line(color = 'NA'),
        panel.grid.minor  = element_line(color ='NA'),
        panel.border = element_rect(fill = 'NA', color = 'black', size = 0.8, linetype = 'solid'),
        axis.text.x =element_text(size = 16, angle = 0,vjust = 0.6,color = 'black',hjust = 0.5,
                                  face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, angle = 0,vjust = 0.6,color = 'black',hjust = 1, face = 'bold'),
        axis.title.y = element_text(vjust=0.2,size = 13))+   theme(legend.position="none", strip.text = element_text(color = "white", size = 14))
micro_p1_all
ggsave("GLM_microbiome_sup.pdf",micro_p1_all,device = "pdf",width = 20,height = 15)

micro_p2_all <-plot_model(micro_glm_model_6, type = "std2", colors = c("#c15560","#7fa4d1"),
                          title = "GLM estimates of effect on microbiome",show.intercept = F,show.legend = T,
                          show.p = TRUE,  grid.breaks=NULL,
                          show.values = T, value.offset = 0.35)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_light()+
  theme(panel.grid = element_line(color = 'NA'),
        panel.grid.minor  = element_line(color ='NA'),
        panel.border = element_rect(fill = 'NA', color = 'black', size = 0.8, linetype = 'solid'),
        axis.text.x =element_text(size = 16, angle = 0,vjust = 0.6,color = 'black',hjust = 0.5,
                                  face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, angle = 0,vjust = 0.25,color = 'black',hjust = 1, face = 'bold'),
        axis.title.y = element_text(vjust=0.2,size = 13))+ theme(legend.position="left")
micro_p2_all 
ggsave("GLM_microbiome.pdf",micro_p2_all,device = "pdf",width = 12,height = 15)

summary(micro_glm_model_6)
micro_summary_text <- capture.output(summary(micro_glm_model_6))
write.csv(micro_summary_text,file = "GLM_microbiome_parameter.csv")
