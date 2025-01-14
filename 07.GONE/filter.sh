for i in TH CH HZ CS YRE; do 
plink --keep $i.pop.txt \
--vcf snps_qc.vcf.gz \
--mac 1 --recode vcf-iid bgz \
--mind 0.6 --chr-set 28 \
--chr 1-28 --allow-extra-chr \
--out $i --thin-count 3000000
done
