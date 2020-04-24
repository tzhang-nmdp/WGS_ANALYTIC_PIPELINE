#!/bin/sh
# 1. basic parameters
outdir=$1
sample_list=$2
con=$3

# 2. run condition in batch
if [ -n "${con}" ] && [ "${con}" == "snpeff" ]; then
        # 2.1 enumerate the sample vcf list
        No2=($(wc -l ${sample_list}))
        for i in $(seq 1 $No2)
        do
                sample=$(awk -v a="$i" '(FNR==a){print $0}' ${sample_list})
                echo ${sample}      
                # 2.2 run snpEff annotation with job number control ( 10 jobs in parallel)
                if (( i%10==0 )); then
                      sh 2.VCF_summary.sh ${sample} ${outdir}
                      else
                      sh 2.VCF_summary.sh ${sample} ${outdir} &                      
                fi
        done
fi
