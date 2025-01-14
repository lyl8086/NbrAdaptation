library(topGO)
library(clusterProfiler)
library(pathview)
library(readr)

setwd("E:/Neosalanx brevirostris/0_R")

term <- read.table("dw_kegg.txt",sep="\t",header=FALSE)
back <- as.character(term$V2)

geneid<- read.csv("HZ_ko.txt",sep="\t",header=FALSE)
id <- as.character(geneid$V2)

result <- enrichKEGG(gene = id, universe = back, keyType = "kegg", organism = 'ko', pvalueCutoff = 1, pAdjustMethod = "none", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
write.csv(summary(result),"KEGG-enrich.csv",row.names =FALSE)


#data(geneList, package="DOSE")

file_name <- paste("4pop_ko","pdf",sep = ".")
pdf(file = file_name)
barplot(result, showCategory = 30)
dev.off()