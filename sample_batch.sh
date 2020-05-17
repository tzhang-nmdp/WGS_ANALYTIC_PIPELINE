#!/bin/sh
# 1. basic parameters
outdir=$1
sample_list=$2
scripts=$3
con=$4
HOME=$(awk '($1~"^HOME"){print $0}' cfg | cut -d "=" -f 2)
SOFTWARE=$(awk '($1~"^SOFTWARE"){print $0}' cfg | cut -d "=" -f 2)
REF=$(awk '($1~"^REF"){print $0}' cfg | cut -d "=" -f 2)

# 2. run annotation condition in batch
if [ -n "${con}" ] && [ "${con}" == "snpeff" ]; then
        # 2.1 enumerate the sample vcf list
        No2=($(wc -l ${outdir}/${sample_list}))
        for i in $(seq 1 $No2)
        do
                sample=$(awk -v a="$i" '(FNR==a){print $0}' ${outdir}/${sample_list})
                echo ${sample} ${con}      
                # 2.2 run snpEff annotation with job number control ( 10 jobs in parallel)
                if (( i%10==0 )); then
                      sh ${HOME}/${scripts} ${outdir} ${sample}
                      else
                      sh ${HOME}/${scripts} ${outdir} ${sample} &                      
                fi
        done
fi

# 2. run vcf combine condition in batch
if [ -n "${con}" ] && [ "${con}" == "vcf_combine" ]; then
        # 2.1 enumerate the sample vcf list
        No2=($(wc -l ${outdir}/${sample_list}))
        rm ${outdir}/${sample_list}.snpeff
        for i in $(seq 1 $No2)
        do
                sample=$(awk -v a="$i" '(FNR==a){print $0}' ${outdir}/${sample_list})
#                echo ${sample} ${con}      
                # 2.2 run snpEff annotation with job number control ( 10 jobs in parallel)
                if (( i%10==0 )); then
                      sh ${HOME}/${scripts} ${outdir} ${sample} ${sample_list}
                      else
                      sh ${HOME}/${scripts} ${outdir} ${sample} ${sample_list} &                      
                fi
        done
     
        # vcf merge by bcftools
        bcftools merge --file-list ${outdir}/${sample_list}.snpeff \
        -f .,PASS --force-samples \
        -m both -Oz \
        --output ${outdir}/${sample_list}.snpeff.combine.vcf.gz \
        --threads 36 
        gunzip -f  ${outdir}/${sample_list}.snpeff.combine.vcf.gz 
fi
