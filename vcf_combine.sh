#!/bin/sh
outdir=$1
sample=$2

(grep ^"#" ${outdir}/${sample}.snpEff.vcf ; grep -v ^"#" ${outdir}/${sample}.snpEff.vcf | sort -T ${HOME}/tmp -k1,1 -k2,2n) | bgzip -c -f > ${outdir}/${sample}.snpEff.vcf.gz
tabix -p vcf -f ${outdir}/${sample}.snpEff.vcf.gz
ls ${outdir}/${sample}.snpEff.vcf.gz >> ${outdir}/${sample_list}.snpeff
