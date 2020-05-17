#!/bin/sh
outdir=$1
sample=$2
sample_list=$3
HOME=$(awk '($1~"^HOME"){print $0}' cfg | cut -d "=" -f 2)
SOFTWARE=$(awk '($1~"^SOFTWARE"){print $0}' cfg | cut -d "=" -f 2)
REF=$(awk '($1~"^REF"){print $0}' cfg | cut -d "=" -f 2)

# sort  and compress the vcf by coordinates
echo ${sample}
(grep ^"#" ${outdir}/snpEff_vcf/${sample}.snpEff.vcf ; \
grep -v ^"#" ${outdir}/snpEff_vcf/${sample}.snpEff.vcf \
| sort -T ${HOME}/tmp -k1,1 -k2,2n) | bgzip -c -f \
> ${outdir}/snpEff_vcf/${sample}.snpEff.vcf.gz

# index the vcf by coordinates
tabix -p vcf -f ${outdir}/snpEff_vcf/${sample}.snpEff.vcf.gz

# make vcf file list
ls ${outdir}/snpEff_vcf/${sample}.snpEff.vcf.gz >> ${outdir}/${sample_list}.snpeff
