library(topGO)

# prepare data.
geneID2GO <- readMappings("dw_go")
geneNames <- names(geneID2GO)
myInterestingGenes <- read.table("fore", header=F, stringsAsFactors=F)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes$V1))
names(geneList) <- geneNames


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
annot = annFUN.gene2GO, gene2GO = geneID2GO)
allGO <- usedGO(object = GOdata) 

# do enrichment.
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher_w <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
# resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, topNodes = length(allGO), numChar=1000, 
	classicFisher = resultFisher,
	weightFisher = resultFisher_w,
	#classicKS = resultKS, 
	#elimKS = resultKS.elim,
	orderBy = "weightFisher", ranksOf = "classicFisher")
fdr <- round(p.adjust(allRes$weightFisher, method="BH"), digits = 5)
allRes <- cbind(allRes, fdr)
# write the result.
write.table(allRes, file="GO_enrichment_Fisher.BP.tsv", quote=F, row.names=F, col.names=T, sep="\t")

# do enrichment.
# resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
# resultFisher_w <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, topNodes = length(allGO), numChar=1000, 
	#classicFisher = resultFisher,
	#weightFisher = resultFisher_w,
	classicKS = resultKS, 
	elimKS = resultKS.elim,
	orderBy = "elimKS", ranksOf = "classicKS")
fdr <- round(p.adjust(allRes$elimKS, method="BH"), digits = 5)
allRes <- cbind(allRes, fdr)
# write the result.
write.table(allRes, file="GO_enrichment_KS.BP.tsv", quote=F, row.names=F, col.names=T, sep="\t")


GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
annot = annFUN.gene2GO, gene2GO = geneID2GO)
allGO <- usedGO(object = GOdata) 

# do enrichment.
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher_w <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
# resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, topNodes = length(allGO), numChar=1000, 
	classicFisher = resultFisher,
	weightFisher = resultFisher_w,
	#classicKS = resultKS, 
	#elimKS = resultKS.elim,
	orderBy = "weightFisher", ranksOf = "classicFisher")
fdr <- round(p.adjust(allRes$weightFisher, method="BH"), digits = 5)
allRes <- cbind(allRes, fdr)
# write the result.
write.table(allRes, file="GO_enrichment_Fisher.CC.tsv", quote=F, row.names=F, col.names=T, sep="\t")

# do enrichment.
# resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
# resultFisher_w <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, topNodes = length(allGO), numChar=1000, 
	#classicFisher = resultFisher,
	#weightFisher = resultFisher_w,
	classicKS = resultKS, 
	elimKS = resultKS.elim,
	orderBy = "elimKS", ranksOf = "classicKS")
fdr <- round(p.adjust(allRes$elimKS, method="BH"), digits = 5)
allRes <- cbind(allRes, fdr)
# write the result.
write.table(allRes, file="GO_enrichment_KS.CC.tsv", quote=F, row.names=F, col.names=T, sep="\t")


GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
annot = annFUN.gene2GO, gene2GO = geneID2GO)
allGO <- usedGO(object = GOdata) 

# do enrichment.
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher_w <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
# resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, topNodes = length(allGO), numChar=1000, 
	classicFisher = resultFisher,
	weightFisher = resultFisher_w,
	#classicKS = resultKS, 
	#elimKS = resultKS.elim,
	orderBy = "weightFisher", ranksOf = "classicFisher")
fdr <- round(p.adjust(allRes$weightFisher, method="BH"), digits = 5)
allRes <- cbind(allRes, fdr)
# write the result.
write.table(allRes, file="GO_enrichment_Fisher.MF.tsv", quote=F, row.names=F, col.names=T, sep="\t")

# do enrichment.
# resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
# resultFisher_w <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, topNodes = length(allGO), numChar=1000, 
	#classicFisher = resultFisher,
	#weightFisher = resultFisher_w,
	classicKS = resultKS, 
	elimKS = resultKS.elim,
	orderBy = "elimKS", ranksOf = "classicKS")
fdr <- round(p.adjust(allRes$elimKS, method="BH"), digits = 5)
allRes <- cbind(allRes, fdr)
# write the result.
write.table(allRes, file="GO_enrichment_KS.MF.tsv", quote=F, row.names=F, col.names=T, sep="\t")