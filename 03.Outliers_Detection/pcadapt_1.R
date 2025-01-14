# PCAdapt

plink --vcf ./HZ_1V2/HZ-1V2-QC.recode.vcf.gz --make-bed --chr-set 28 --allow-extra-chr --out test

library(pcadapt)
library(qvalue)
filename <- read.pcadapt('./Pcadapt/HZ-1V2.bed', type = "bed")
res <- pcadapt(filename, K = 20, LD.clumping = list(size = 50, thr = 0.1))
pdf("screeplot_rm.pdf")
plot(res, option = "screeplot")
dev.off()
k<-2
res <- pcadapt(filename, K = as.numeric(k), LD.clumping = list(size = 50, thr = 0.1))
qval <- qvalue(res$pvalues)$qvalues
write.table(res$pvalues, "pcadapt.allp", col.names=F, row.names=F, quote=F)
outliers <- which(qval <= 0.1)
write.table(outliers, "pcadapt.outliers.id", col.names=F, row.names=F, quote=F)
