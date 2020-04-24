#!/bin/sh
# 1. basic parameters
outdir=$1
sample_list=$2
snpeff_script=$3
con=$4

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
                      sh ${snpeff_script} ${outdir} ${sample}
                      else
                      sh ${snpeff_script} ${outdir} ${sample} &                      
                fi
        done
fi
