
setwd("D:/Desktop/Projects/sewage/Figure/Final_Figure/Figure5/IS/sharing/sharing_real/")

library(ggplot2)
library(gggenes)

a <- read.delim("genome_combined_AAC(6')-Ib7.txt", sep = "\t")
dummies_a <- make_alignment_dummies(a,  
                                  aes(xmin = start, xmax = end, y = molecule, id = gene),
                                  on = "mphE")
p_a <- ggplot(a, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
  geom_gene_arrow() +
  theme_genes()
p_a
ggsave("AAC(6')-Ib7.pdf",p_a,device = "pdf",width = 7.5,height = 5)

b <- read.delim("genome_combined_AAC(6')-Ib-cr11.txt", sep = "\t")
dummies_b <- make_alignment_dummies(b,  
                                    aes(xmin = start, xmax = end, y = molecule, id = gene),
                                    on = "mphE")
p_b <- ggplot(b, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
  geom_gene_arrow() +
  theme_genes()
p_b
ggsave("AAC(6')-Ib-cr11.pdf",p_b,device = "pdf",width = 7.5,height = 5)

c <- read.delim("genome_combined_mphE.txt", sep = "\t")
dummies_c <- make_alignment_dummies(c,  
                                    aes(xmin = start, xmax = end, y = molecule, id = gene),
                                    on = "mphE")
p_c <- ggplot(c, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
  geom_gene_arrow() +
  theme_genes()
p_c
ggsave("mphE.pdf",p_c,device = "pdf",width = 7.5,height = 5)

d <- read.delim("genome_combined_sul1.txt", sep = "\t")
dummies_d <- make_alignment_dummies(d,  
                                    aes(xmin = start, xmax = end, y = molecule, id = gene),
                                    on = "mphE")
p_d <- ggplot(d, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
  geom_gene_arrow() +
  theme_genes()
p_d
ggsave("sul1.pdf",p_d,device = "pdf",width = 7.5,height = 5)

e <- read.delim("genome_combined_sul2.txt", sep = "\t")
dummies_e <- make_alignment_dummies(e,  
                                    aes(xmin = start, xmax = end, y = molecule, id = gene),
                                    on = "mphE")
p_e <- ggplot(e, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
  geom_gene_arrow() +
  theme_genes()
p_e
ggsave("sul2.pdf",p_e,device = "pdf",width = 7.5,height = 5)

f <- read.delim("genome_combined_tet(X3).txt", sep = "\t")
dummies_f <- make_alignment_dummies(f,  
                                    aes(xmin = start, xmax = end, y = molecule, id = gene),
                                    on = "mphE")
p_f <- ggplot(f, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
  geom_gene_arrow() +
  theme_genes()
p_f
ggsave("tet(X3).pdf",p_f,device = "pdf",width = 7.5,height = 5)
