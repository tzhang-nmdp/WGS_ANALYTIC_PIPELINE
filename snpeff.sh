#!/bin/sh
# 1. parameter setting #
HOME=$(awk '($1~"^HOME"){print $0}' cfg | cut -d "=" -f 2)
SOFTWARE=$(awk '($1~"^SOFTWARE"){print $0}' cfg | cut -d "=" -f 2)
REF=$(awk '($1~"^REF"){print $0}' cfg | cut -d "=" -f 2)
outdir=$1
sample=$2
mkdir ${outdir}/snpEff_vcf

# 2. snpEff basic variant annotation
echo -e "${vcf} snpEff startTime::::::::\c"; date

java -d64 -Xmx3g -jar ${SOFTWARE}/snpEff/snpEff.jar \
-v GRCh38.86 ${outdir}/${sample}.vcf \
-canon -stats ${outdir}/snpEff_vcf/${sample_vcf}.html \
> ${outdir}/snpEff_vcf/${sample}.snpEff.vcf

echo -e "${vcf} snpEff endTime::::::::\c"; date
