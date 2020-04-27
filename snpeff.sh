#!/bin/sh
# 1. parameter setting #
HOME=/home/tzhang
outdir=$1
sample_vcf=$2
mkdir ${outdir}/snpEff_vcf

# 2. snpEff basic variant annotation
echo -e "${vcf} snpEff startTime::::::::\c"; date

java -d64 -Xmx3g -jar ${HOME}/software/snpEff/snpEff.jar \
-v GRCh38.86 ${outdir}/snpEff_vcf/${sample_vcf} \
-canon -stats ${outdir}/snpEff_vcf/${sample_vcf}.html \
> ${outdir}/snpEff_vcf/${sample_vcf}.snpEff.vcf

echo -e "${vcf} snpEff endTime::::::::\c"; date
