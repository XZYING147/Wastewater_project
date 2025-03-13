rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure4/GAM")

library(mgcv)

virus_data <- read.csv("GLM_sewage.csv",check.names = F, header = T, row.names = 1)

resistome_gam_model <- gam(Resistome_Shannon ~ s(Microbiome_Shannon) + s(Color_Fold) + 
                     s(pH) + s(Total_Suspended_Solid),
                   data = virus_data, family = gaussian, select=TRUE)

summary(resistome_gam_model)

svg("GAM_resistome_microbiome_shannon.svg", width = 7.5, height = 5)
plot(resistome_gam_model, select = 1, pch = 20, shade = TRUE, residuals = TRUE)
dev.off() 

svg("GAM_resistome_colorfold.svg", width = 7.5, height = 5)
plot(resistome_gam_model, select = 2, pch = 20, shade = TRUE, residuals = TRUE)
dev.off() 

svg("GAM_resistome_ph.svg", width = 7.5, height = 5)
plot(resistome_gam_model, select = 3, pch = 20, shade = TRUE, residuals = TRUE)
dev.off() 

svg("GAM_resistome_TSS.svg", width = 7.5, height = 5)
plot(resistome_gam_model, select = 4, pch = 20, shade = TRUE, residuals = TRUE)
dev.off() 


micro_gam_model <- gam(Microbiome_Shannon ~ s(Resistome_Shannon) + s(pH),
                         data = virus_data, family = gaussian, select=TRUE)

summary(micro_gam_model)

svg("GAM_microbiome_resistome_shannon.svg", width = 7.5, height = 5)
plot(micro_gam_model, select = 1, pch = 20, shade = TRUE, residuals = TRUE)
dev.off() 

svg("GAM_microbiome_ph.svg", width = 7.5, height = 5)
plot(micro_gam_model, select = 2, pch = 20, shade = TRUE, residuals = TRUE)
dev.off() 
