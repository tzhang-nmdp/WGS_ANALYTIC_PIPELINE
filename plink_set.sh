#!/bin/sh
# start with snpEff annotated vcf
outdir=$1
vcf=$2
chr=$3
temp=$4
echo ${chr}
HOME=/home/tzhang

echo "#1. chromosome variant extraction" && date
awk -v a="$chr" '($1~"#")||($1==a){print $0}' ${outdir}/${vcf}  > ${outdir}/${vcf}.${chr}

echo "#2 chromosome variant transformation" && date
if [ ${temp} == "somatic" ] ; then 
echo "somatic"
python ${HOME}/vcf_id.py ${outdir}/${vcf}.${chr} somatic > ${outdir}/${vcf}.${chr}.vcf
elif [ ${temp} == "germline" ] ; then 
echo "germline"
python ${HOME}/vcf_id.py ${outdir}/${vcf}.${chr} germ > ${outdir}/${vcf}.${chr}.vcf
fi

# check vcf transformation
if [ $? -eq 0 ]; then
echo "vcf file ready"
rm ${outdir}/${vcf}.${chr}
else
echo ${vcf} >> ${outdir}/vcf.error
exit
fi

echo "#3. chromosome gene extraction" && date
python ${HOME}/plink_set.py \
${outdir}/${vcf}.${chr}.vcf \
${outdir}/${vcf}.${chr}.vcf.geneset \
${HOME}/database/database_port/gnomad_ab0.05_${chr} 

cat ${outdir}/${vcf}.${chr}.vcf.geneset[._]* > ${outdir}/${vcf}.${chr}.vcf.geneset
awk 'FNR>1{split($8,b,";");split(b[1],a,"="); if ($1!~"#") print a[2]"\t"$3}' ${outdir}/${vcf}.${chr}.vcf >  ${outdir}/${vcf}.${chr}.SKAT

echo "#4. for SKAT geneset"
awk 'FNR>1{split($8,b,";");split(b[1],a,"="); if ($1!~"#") print a[2]"\t"$3}' ${outdir}/${vcf}.${chr}.vcf.gene >  ${outdir}/${vcf}.${chr}.SKAT
