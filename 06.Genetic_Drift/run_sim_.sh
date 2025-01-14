#!
g=$1
af=$2
rep=1000
mkdir ${g}Ge
cp $af genetic_drift_fast.py ${g}Ge/

pushd ${g}Ge
echo "Nyr	2451	$g
Nth	4919	$g
Nhz	1600	$g
Nch	4021	$g
Ncs	327	$g" >Ne.txt

# generate cmd
while read line; do 
    ne=`echo $line | awk '{print $2}'`;
    g=`echo $line | awk '{print $3}'`;
    n=`echo $line | awk '{print $1}'`;
    for i in `seq 1 $rep`; do
        echo "python3 genetic_drift_fast.py $af $ne $g | gzip > $n.$i.txt.gz" >>run.cmd;
    done
done < Ne.txt
# run！
ParaFly -c run.cmd -CPU 128

# parse results
mkdir change
for p in Nth Nhz Nch Ncs; 
do for i in `seq 1 $rep`; 
	do paste <(zcat Nyr.$i.txt.gz) <(zcat $p.$i.txt.gz) | awk -v OFS="\t" '{print $1,$2-$1}' | gzip >change/$p.$i.txt.gz;
	done; 
done

pushd change
for i in `seq 1 $rep`; do paste <(zcat Nth.$i.txt.gz ) <(zcat Nhz.$i.txt.gz | cut -f2) <(zcat Nch.$i.txt.gz | cut -f2) <(zcat Ncs.$i.txt.gz | cut -f2) >$i.change.txt; done

rm -rf Nth.*.txt.gz Nhz.*.txt.gz Nch.*.txt.gz Ncs.*.txt.gz

cat << EOF >plot.R
library(ggpubr)
library(reshape2)
library(data.table)

plist <- list()
for (i in 1:$rep) {
    x<-fread(paste0(i,".change.txt"))
    x<-x[,-1]
    names(x)<-c("TH","HZ","CH","CS")
    y<-melt(x)
    names(y)<-c("Pairs","DeltaF")
    nn<-''
    for (j in seq(0,0.7,0.05)) {
        m<-nrow(x[TH>j & HZ>j & CH>j & CS>j,])
        m<-m+nrow(x[TH<(-j) & HZ<(-j) & CH<(-j) & CS<(-j),])
        n<-paste0(j,": ",m)
        nn<<-paste(n,nn,sep="\n")
    }
    
    plist[[i]]<-gghistogram(y, x = "DeltaF", color = "Pairs", fill = "Pairs", alpha=0.5) +
    annotate(geom="text", x=0.3, y=5000, label=nn, size=2.5,
              color="red")
}

pdf("hist.pdf", height=40,width=20)
ggarrange(plotlist = plist, ncol=5, nrow=20,common.legend = TRUE)
EOF

Rscript plot.R

