
vcf=$1
rep=20
out=${vcf%%.vcf.gz}
plink --vcf $vcf --recode 12 --make-bed --out $out --allow-extra-chr --autosome-num 77

bioawk -t \
'{$2=NR;
sub(/^I$/,"1",$1);
sub(/^II$/,"2",$1);
sub(/^III$/,"3",$1);
sub(/^IV$/,"4",$1);
sub(/^IX$/,"9",$1);
sub(/^V$/,"5",$1);
sub(/^VI$/,"6",$1);
sub(/^VII$/,"7",$1);
sub(/^VIII$/,"8",$1);
sub(/^X$/,"10",$1);
sub(/^XI$/,"11",$1);
sub(/^XII$/,"12",$1);
sub(/^XIII$/,"13",$1);
sub(/^XIV$/,"14",$1);
sub(/^XIX$/,"19",$1);
sub(/^XV$/,"15",$1);
sub(/^XVI$/,"16",$1);
sub(/^XVII$/,"17",$1);
sub(/^XVIII$/,"18",$1);
sub(/^XX$/,"20",$1);
sub(/^XXI$/,"21",$1);
sub(/^XXII$/,"22",$1);
sub(/^XXIII$/,"23",$1);
sub(/^XXIV$/,"24",$1);
print}' \
$out.map >map.1
mv map.1 $out.map

:>run.cmd
for x in `seq 1 $rep`; do
    cp -a /opt/bio/populations/GONE/Linux GONE_$x
    cp $out.ped $out.map GONE_$x
    echo "pushd GONE_$x; bash script_GONE.sh $out" >>run.cmd
done

echo "
library(ggpubr)
library(data.table)

x<-fread(paste0(\"GONE_1\",\"/Output_Ne_$out\"), skip=2, header=F)
names(x)<-c(\"Gen\",\"Ne\")
p<-ggplot(data=x[1:200,], aes(x=Gen, y=Ne))+
geom_line(alpha=0.2)+theme_pubr()+
labs(title = \"$out, cMMb=1, hc=0.05\")+
theme(plot.title=element_text(color=\"red\",hjust=0.5))
a<-fread(paste0(\"GONE_1\",\"/Output_Ne_$out\"), skip=2, header=F)
names(a)<-c(\"Gen\",\"rep1\")

for (i in seq(2,$rep)) {
	x<-fread(paste0(\"GONE_\",i,\"/Output_Ne_$out\"), skip=2, header=F)
	names(x)<-c(\"Gen\",\"Ne\")
    a<<-a[x,on=\"Gen\"]
    names(a)[i+1]<-paste0(\"rep\",i)
	p<<-p+geom_line(data=x[1:200,],aes(x=Gen, y=Ne),alpha=0.2)
}

mm <- function(x) { exp(mean(log(x)))}
a[,Ne := apply(a[,2:($rep+1)],1,mm)]
fwrite(a,\"$out.xls\",sep=\"\t\",quote=FALSE)
p<-p+geom_line(data=a[1:200,],aes(x=Gen, y=Ne),color=\"red\",alpha=0.5, size=1)
ggsave(\"$out.GONE.pdf\",p)
" >$out.plot.R
