#!/bin/sh
############################################################################################################################################################
# example commandline sh sample_batch.sh -s sample.list -o example -c set -s
############################################################################################################################################################

#------------------------------------------------------------------------------------------------------------------------------
## Default parameter option
#------------------------------------------------------------------------------------------------------------------------------

HOME=$(awk '($1~"^HOME"){print $0}' cfg | cut -d "=" -f 2)
SOFTWARE=$(awk '($1~"^SOFTWARE"){print $0}' cfg | cut -d "=" -f 2)
REF=$(awk '($1~"^REF"){print $0}' cfg | cut -d "=" -f 2)

#------------------------------------------------------------------------------------------------------------------------------
## Option setup for paramters 
#------------------------------------------------------------------------------------------------------------------------------
while getopts ":s:o:c:s:" opt
do
    case $opt in
        s)
            sample_list=$OPTARG
        ;;
        o)
            outdir=$OPTARG
        ;;
        c)
            con=$OPTARG
        ;;
        ?)
            echo "unknown argument"
            exit 1
        ;;
    esac
done
if [ ! ${sample_vcf} ]; then
    echo "-s invalid"
    exit 1
fi

if [ ! ${oPath} ]; then
    echo "-o invalid"
    exit 1
fi

if [ ! ${con} ]; then
    echo "-c invalid"
    exit 1
fi

date=$(date +%F)
shellPath=$(cd "$(dirname "$0")"; pwd)
echo ${shellPath}
echo ${date}

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
