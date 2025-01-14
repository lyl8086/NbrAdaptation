rep=$1

for i in `seq 1 $rep`; do
    
    for j in 0.1 0.15 0.2 0.25; do
        
        cat << EOF >>run.cmd
awk -v a=$j -v b=-$j '(\$2>a && \$3>a && \$4>a && \$5>a) || (\$2<b && \$3<b && \$4<b && \$5<b) {c++} END {print c}' $i.change.txt >>$j.txt
EOF
     done
done

ParaFly -c run.cmd -CPU 128

cat << EOF >change.hist.R
library(ggpubr)
library(data.table)

plist <- list()
for (i in 0.1 0.15 0.2 0.25) {
    x<-fread(paste0(i,".txt"))
    plist[[i]]<-gghistogram(x, x = "V1", color = "gray", fill = "gray", xlab="DeltaF")
}

pdf("change_hist.pdf", height=7,width=7)
ggarrange(plotlist = plist, ncol=2, nrow=2,common.legend = TRUE)
EOF

Rscript change.hist.R