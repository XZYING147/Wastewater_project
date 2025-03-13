rm(list = ls())
setwd("D:/Desktop/Projects/wastewater/Figure/Figure2/Upset_network")

library(UpSetR)
outFile="intersectGenes.txt"
outPic="Upset_network_resistome.pdf"

files=dir()
files=grep("txt$",files,value=T)
geneList=list()

for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F)
  geneNames=as.vector(rt[,1])
  geneNames=gsub("^ | $","",geneNames)
  uniqGene=unique(geneNames)
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}

#绘制UpSetͼ
upsetData=fromList(geneList)
pdf(file=outPic,onefile = FALSE,width=9,height=6)
upset(upsetData,
      nsets = length(geneList),
      nintersects = 50,
      order.by = "freq",
      show.numbers = "yes",
      point.size = 2,
      matrix.color="black",
      line.size = 0.8,
      mainbar.y.label = "Intersections", 
      sets.x.label = "Number")
dev.off()