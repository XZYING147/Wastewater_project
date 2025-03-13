rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure2/Inner_cor")
options(stringsAsFactors=F)

library(corrplot)
library(dplyr)
library(reshape2)
library(igraph)
library(ggraph)

rt <- read.delim("rpkm.type_sample.txt",row.names = 1)
rt_t <- t(rt)
rt_final <- as.data.frame(rt_t)
rt_fli <- rt_final[1:22,1:22]
rt_hos <- rt_final[23:38,1:22]
rt_res <- rt_final[39:93,1:22]
rt_wtp <- rt_final[94:252,1:22]

rt_final <- rt_final[ , !(names(rt_final) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]
rt_fli <- rt_fli[ , !(names(rt_fli) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]
rt_res <- rt_res[ , !(names(rt_res) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]
rt_hos <- rt_hos[ , !(names(rt_hos) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]
rt_wtp <- rt_wtp[ , !(names(rt_wtp) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]

cor1=cor(log10(rt_final+1))
similarity_long <- melt(cor1)
similarity_long_filtered <- similarity_long %>% filter(value != 1)
similarity_long_filtered <- similarity_long_filtered %>% filter(abs(value) > 0.5)
similarity_long_filtered$width <- abs(similarity_long_filtered$value)*1.5


nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges <- similarity_long_filtered

net<-graph_from_data_frame(d=edges,vertices = nodes,directed = T)

p_all <- ggraph(net, layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ 
  # geom_node_text(aes(label=Name))+
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) + 
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right" 
  )
ggsave("cor_inner_all.pdf",p_all,device = "pdf",width = 7.5,height = 6.5)


cor_fli=cor(log10(rt_fli+1))
similarity_long_fli <- melt(cor_fli)
similarity_long_filtered_fli <- similarity_long_fli %>% filter(value != 1)
similarity_long_filtered_fli <- similarity_long_filtered_fli %>% filter(abs(value) > 0.5)
similarity_long_filtered_fli$width <- abs(similarity_long_filtered_fli$value)*1.5

nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges_fli <- similarity_long_filtered_fli

net_fli<-graph_from_data_frame(d=edges_fli,vertices = nodes,directed = T)

p_fli <- ggraph(net_fli, layout = "linear", circular = TRUE)+
  # geom_edge_link(aes(start_cap = label_rect(node1.name),
  #                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value)
  # )+ #arrow = arrow(length = unit(2, 'mm'))
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ #,color=Group
  # geom_node_text(aes(label=Name))+
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) + 
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right"
  )
ggsave("cor_inner_fli.pdf",p_fli,device = "pdf",width = 7.5,height = 6.5)


cor_res=cor(log10(rt_res+1))
similarity_long_res <- melt(cor_res)
similarity_long_filtered_res <- similarity_long_res %>% filter(value != 1)
similarity_long_filtered_res <- similarity_long_filtered_res %>% filter(abs(value) > 0.5)
similarity_long_filtered_res$width <- abs(similarity_long_filtered_res$value)*1.5

nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges_res <- similarity_long_filtered_res

net_res<-graph_from_data_frame(d=edges_res,vertices = nodes,directed = T)

p_res <- ggraph(net_res, layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ 
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) + 
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right"
  )
ggsave("cor_inner_res.pdf",p_res,device = "pdf",width = 7.5,height = 6.5)


cor_hos=cor(log10(rt_hos+1))
similarity_long_hos <- melt(cor_hos)
similarity_long_filtered_hos <- similarity_long_hos %>% filter(value != 1)
similarity_long_filtered_hos <- similarity_long_filtered_hos %>% filter(abs(value) > 0.5)
similarity_long_filtered_hos$width <- abs(similarity_long_filtered_hos$value)*1.5

nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges_hos <- similarity_long_filtered_hos

net_hos<-graph_from_data_frame(d=edges_hos,vertices = nodes,directed = T)

p_hos <- ggraph(net_hos, layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ 
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) + 
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right"
  )
ggsave("cor_inner_hos.pdf",p_hos,device = "pdf",width = 7.5,height = 6.5)


cor_wtp=cor(log10(rt_wtp+1))
similarity_long_wtp <- melt(cor_wtp)
similarity_long_filtered_wtp <- similarity_long_wtp %>% filter(value != 1)
similarity_long_filtered_wtp <- similarity_long_filtered_wtp %>% filter(abs(value) > 0.5)
similarity_long_filtered_wtp$width <- abs(similarity_long_filtered_wtp$value)*1.5

nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges_wtp <- similarity_long_filtered_wtp

net_wtp<-graph_from_data_frame(d=edges_wtp,vertices = nodes,directed = T)

p_wtp <- ggraph(net_wtp, layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ 
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) +
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right"
  )
ggsave("cor_inner_wtp.pdf",p_wtp,device = "pdf",width = 7.5,height = 6.5)


rt_time <- read.delim("rpkm.type_sample_time.txt",row.names = 1)
rt_time_t <- t(rt_time)
rt_time_final <- as.data.frame(rt_time_t)
rt_feb <- rt_time_final[1:40,1:22]
rt_mar <- rt_time_final[41:193,1:22]
rt_apr <- rt_time_final[194:252,1:22]

rt_feb <- rt_feb[ , !(names(rt_feb) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]
rt_mar <- rt_mar[ , !(names(rt_mar) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]
rt_apr <- rt_apr[ , !(names(rt_apr) %in% c("tunicamycin","tetracenomycin_C","puromycin","other_peptide_antibiotics","fusidic_acid","defensin","bicyclomycin","multidrug"))]

cor_feb=cor(log10(rt_feb+1))
similarity_long_feb <- melt(cor_feb)
similarity_long_filtered_feb <- similarity_long_feb %>% filter(value != 1)
similarity_long_filtered_feb <- similarity_long_filtered_feb %>% filter(abs(value) > 0.5)
similarity_long_filtered_feb$width <- abs(similarity_long_filtered_feb$value)*1.5


nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges_feb <- similarity_long_filtered_feb

net_feb<-graph_from_data_frame(d=edges_feb,vertices = nodes,directed = T)

p_feb <- ggraph(net_feb, layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ 
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) + 
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right"
  )
ggsave("cor_inner_feb.pdf",p_feb,device = "pdf",width = 7.5,height = 6.5)


cor_mar=cor(log10(rt_mar+1))
similarity_long_mar <- melt(cor_mar)
similarity_long_filtered_mar <- similarity_long_mar %>% filter(value != 1)
similarity_long_filtered_mar <- similarity_long_filtered_mar %>% filter(abs(value) > 0.5)
similarity_long_filtered_mar$width <- abs(similarity_long_filtered_mar$value)*1.5

nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges_mar <- similarity_long_filtered_mar

net_mar<-graph_from_data_frame(d=edges_mar,vertices = nodes,directed = T)

p_mar <- ggraph(net_mar, layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ 
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) + 
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right"
  )
ggsave("cor_inner_mar.pdf",p_mar,device = "pdf",width = 7.5,height = 6.5)


cor_apr=cor(log10(rt_apr+1))
similarity_long_apr <- melt(cor_apr)
similarity_long_filtered_apr <- similarity_long_apr %>% filter(value != 1)
similarity_long_filtered_apr <- similarity_long_filtered_apr %>% filter(abs(value) > 0.5)
similarity_long_filtered_apr$width <- abs(similarity_long_filtered_apr$value)*1.5

nodes <- read.table(file = "nodes.txt", sep = "\t", header = T, check.names = FALSE)
edges_apr <- similarity_long_filtered_apr

net_apr<-graph_from_data_frame(d=edges_apr,vertices = nodes,directed = T)

p_apr <- ggraph(net_apr, layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name),edge_width=width, edge_color=value))+ 
  coord_fixed()+
  geom_node_point(aes(size=Size,color=Group),shape = 16,show.legend = FALSE)+ 
  scale_size_area(max_size = 15, name = NULL, guide =FALSE)+
  scale_color_manual(values = c("Aminoglycoside"="#f4c0bd", "Antibacterial_fatty_acid"="#ead1d1", "Bacitracin"="#d1eee9","Beta_lactam"="#f7f4d7","Bleomycin"="#efd094",
                                "Chloramphenicol"="#bfdfd2","Florfenicol"="#e3edae","Fosfomycin"="#ccdef4","MLS"="#f6e2eb","Mupirocin"="#f4ecb4",
                                "Novobiocin"="#c2e9e2","Pleuromutilin_tiamulin"="#f6dbf1","Polymyxin"="#a1c0c7","Quinolone"="#e8ddaf","Rifamycin"="#d2d2d2",
                                "Streptothricin"="#718991","Sulfonamide"="#b8b3d6","Tetracycline"="#a9cbdf","Trimethoprim"="#f5dfe1","Vancomycin"="#aec6da")) + 
  scale_edge_color_gradient2(low = "#5ebcc2",mid="white",high = "#5f5f5d")+ 
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x =element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right"
  )
ggsave("cor_inner_apr.pdf",p_apr,device = "pdf",width = 7.5,height = 6.5)