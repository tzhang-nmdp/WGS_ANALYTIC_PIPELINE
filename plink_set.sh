#!/bin/sh
# start with snpEff annotated vcf
vcf=$1
outdir=$2
chr=$3
con=$4
echo ${chr}

echo "#1. chromosome variant extraction" && date
awk -v a="$chr" '($1~"#")||($1==a){print $0}' ${outdir}/${vcf}  > ${outdir}/${vcf}.${chr}

echo "#2 chromosome variant transformation" && date
if [ ${con} == "somatic" ] ; then 
echo "somatic"
python /home/tzhang/MDS_data/vcf_id.py ${outdir}/${vcf}.${chr} somatic > ${outdir}/${vcf}.${chr}.vcf
elif [ ${con} == "germline" ] ; then 
echo "germline"
python /home/tzhang/MDS_data/vcf_id.py ${outdir}/${vcf}.${chr} germ > ${outdir}/${vcf}.${chr}.vcf
fi

echo "#3. chromosome gene extraction" && date
python /home/tzhang/MDS_data/plink_set.py ${outdir}/${vcf}.${chr}.vcf ${outdir}/${vcf}.${chr}.vcf.geneset /home/tzhang/database/database_port/gnomad_ab0.05_${chr} 
cat ${outdir}/${vcf}.${chr}.vcf.geneset[._]* > ${outdir}/${vcf}.${chr}.vcf.geneset
awk 'FNR>1{split($8,b,";");split(b[1],a,"="); if ($1!~"#") print a[2]"\t"$3}' ${outdir}/${vcf}.${chr}.vcf >  ${outdir}/${vcf}.${chr}.SKAT

echo "#4. for SKAT geneset"
awk 'FNR>1{split($8,b,";");split(b[1],a,"="); if ($1!~"#") print a[2]"\t"$3}' ${outdir}/${vcf}.${chr}.vcf.gene >  ${outdir}/${vcf}.${chr}.SKAT
