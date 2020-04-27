#!/bin/sh
outdir=$1
sample=$2
HOME=/home/tzhang

# sort  and compress the vcf by coordinates
(grep ^"#" ${outdir}/${sample}.snpEff.vcf ; grep -v ^"#" ${outdir}/${sample}.snpEff.vcf \
| sort -T ${HOME}/tmp -k1,1 -k2,2n) \ 
| bgzip -c -f \
> ${outdir}/${sample}.snpEff.vcf.gz

# index the vcf by coordinates
tabix -p vcf -f ${outdir}/${sample}.snpEff.vcf.gz

# make vcf file list
ls ${outdir}/${sample}.snpEff.vcf.gz >> ${outdir}/${sample_list}.snpeff
