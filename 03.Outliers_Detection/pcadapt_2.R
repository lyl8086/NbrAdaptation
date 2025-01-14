tmp_vcf<-readLines("/mnt/DATA/DAT2/home/yanghao/FWA/TH_YRE/TH_YRE.vcf")
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
tmp_vcf_data<-read.table("/mnt/DATA/DAT2/home/yanghao/FWA/TH_YRE/TH_YRE.vcf", stringsAsFactors = FALSE)
names(tmp_vcf_data)<-vcf_names
id<-readLines("/mnt/DATA/DAT2/home/yanghao/FWA/TH_YRE/PCAdapt/TH_YRE.pcadapt.outliers.id")
outliers<-data.frame()
for (i in id){
  outliers<-rbind(outliers,tmp_vcf_data[i,])
}
write.table(outliers, "TH_YRE.outliers.vcf", col.names=F, row.names=F, quote=F,sep = "\t")