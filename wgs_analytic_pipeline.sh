#!/bin/sh
############################################################################################################################################################
# example commandline sh ../../12.plink_gene_set.assoc.sh -s NMDP.somatic.vcf.R -o /home/tzhang/MDS_data/test -c set -t somatic -n 1 22
############################################################################################################################################################

#------------------------------------------------------------------------------------------------------------------------------
## Default parameter option
#------------------------------------------------------------------------------------------------------------------------------

HOME=/home/tzhang
PLINK=${HOME}/software/plink_1.9/plink

#------------------------------------------------------------------------------------------------------------------------------
## Option setup for paramters 
#------------------------------------------------------------------------------------------------------------------------------
while getopts ":s:o:c:t:n:" opt
do
    case $opt in
        s)
            sample_vcf=$OPTARG
        ;;
        o)
            oPath=$OPTARG
        ;;
        c)
            con=$OPTARG
        ;;
        t)
            temp=$OPTARG
        ;;
        n)
            No=$OPTARG
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

if [ ! ${temp} ]; then
    echo "-t invalid"
    exit 1
fi

if [ ! ${No} ]; then
    echo "-n invalid"
    exit 1
fi

No1=$(echo ${No} | awk '{print $1}')
No2=$(echo ${No} | awk '{print $2}')

date=$(date +%F)
shellPath=$(cd "$(dirname "$0")"; pwd)
echo ${shellPath}

#------------------------------------------------------------------------------------------------------------------------------
## vcf data transformation
#------------------------------------------------------------------------------------------------------------------------------

if [ -n "${con}" ] && ([ "${con}" == "prep" ] || [ "${con}" == "all" ]); then
    echo "############################################## vcf data transformation ##################################"
    echo -e "${sample} plink-pre startTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    for x in $(seq $No1 $No2) X
        do
        chr="chr"${x}
        echo ${sample_vcf} " go to " ${chr}
        sh ${HOME}/plink_set.sh ${oPath} ${sample_vcf} ${chr} ${temp}
    done
    echo -e "${sample} plink-pre endTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
fi


#------------------------------------------------------------------------------------------------------------------------------
## variant PLINK association test
#------------------------------------------------------------------------------------------------------------------------------

if [ -n "${con}" ] && ([ "${con}" == "sgl" ]|| [ "${con}" == "all" ]); then
    echo "############################################## variant PLINK association test ##################################"
    echo -e "${sample} plink-single startTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    mkdir  ${oPath}/sgl
    for x in $(seq $No1 $No2) X
    do
        chr="chr"$x
        echo ${sample_vcf} " go to " ${chr}
        ${PLINK} --vcf ${oPath}/${sample_vcf}.${chr}.vcf --threads 128 --pheno ${HOME}/${sample_vcf}.fam --mpheno 4 --allow-no-sex --covar ${HOME}/${sample_vcf}.cov --covar-number 2-13 --logistic perm --pfilter 0.05 --out ${oPath}/sgl/${sample_vcf}.${chr}.vcf.sig.R
        ${PLINK} --vcf ${oPath}/${sample_vcf}.${chr}.vcf --threads 128 --pheno ${HOME}/${sample_vcf}.fam --mpheno 4 --allow-no-sex --covar ${HOME}/${sample_vcf}.cov --covar-number 2-13 --assoc perm --pfilter 0.05 --out ${oPath}/sgl/${sample_vcf}.${chr}.vcf.sig.R
    done
    rm ${oPath}/sgl/*.log
    rm ${oPath}/sgl/*sex
    echo -e "${sample} plink-single endTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
fi

#------------------------------------------------------------------------------------------------------------------------------
## variant conditional logistic regression association test
#------------------------------------------------------------------------------------------------------------------------------

if [ -n "${con}" ] && ([ "${con}" == "clr" ]|| [ "${con}" == "all" ]); then
    echo "############################################## variant conditional logistic regression association test ##################################"
    echo -e "${sample} clr startTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    mkdir  ${oPath}/clr
    for x in $(seq $No1 $No2) X
    do
        chr="chr"$x
        echo ${sample_vcf} " go to " ${chr}
        awk '($1~"#CHROM")||($1~"CHROM"){print $0}' ${oPath}/${sample_vcf}.${chr}.vcf  | sed 's%#%%g' > ${oPath}/${sample_vcf}.header.tmp
        #sed 's%#CHROM%CHROM%g' ${oPath}/${sample_vcf}.${chr}.vcf > ${oPath}/${sample_vcf}.${chr}.vcf.gene.tmp
        Rscript ${HOME}/R_CLOGIT.R -i ${oPath}/${sample_vcf}.${chr}.vcf -s ${HOME}/IDR_pair_outcome.txt -o ${oPath}/clr/${sample_vcf}.${chr}.vcf.R.clr -d ${oPath}/${sample_vcf}.header.tmp
    done
    echo -e "${sample} clr endTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log1
fi

#------------------------------------------------------------------------------------------------------------------------------
## PLINK gene burden test
#------------------------------------------------------------------------------------------------------------------------------

if [ -n "${con}" ] && ([ "${con}" == "set" ]|| [ "${con}" == "all" ]); then
    echo "#############################################PLINK-GENE-SET-ASSOCIATION-TEST-ASSO##################################"
    echo -e "${sample} plink-set startTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    mkdir ${oPath}/plink_gb
    for x in $(seq $No1 $No2) X
    do
        chr="chr"$x
        echo ${sample} " go to " ${chr}
        No=$(ls ${oPath}/${sample_vcf}.${chr}.vcf.geneset_* | wc -l)
        No3=$((No-1))
        
        ${PLINK} --vcf ${oPath}/${sample_vcf}.${chr}.vcf --threads 128 --pheno ${HOME}/${sample_vcf}.fam --mpheno 4 --allow-no-sex --recode --make-bed --out ${oPath}/${sample_vcf}.${chr}.R
        for i in $(seq 0 $No3)
        do
        ${PLINK} --bfile ${oPath}/${sample_vcf}.${chr}.R --pheno ${HOME}/${sample_vcf}.fam --mpheno 4  --covar ${HOME}/${sample_vcf}.cov --covar-number 2-13 --threads 128  --allow-no-sex --set ${oPath}/${sample_vcf}.${chr}.vcf.geneset_${i} --logistic --assoc perm set-test --set-p 0.05  --pfilter 0.05 --out ${oPath}/plink_gb/${sample_vcf}.${chr}.vcf.set.R_${i} 
        done
        echo -e "${sample} plink-set endTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    done
    rm ${oPath}/plink_gb/*.log
    rm ${oPath}/plink_gb/*sex
fi

#------------------------------------------------------------------------------------------------------------------------------
## SKAT-O gene association test
#------------------------------------------------------------------------------------------------------------------------------

if [ -n "${con}" ] && ([ "${con}" == "skat" ]|| [ "${con}" == "all" ]); then
    echo "############################################# SKAT-O gene association test ##################################"
    echo -e "${sample} skat-o startTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    mkdir -p ${oPath}/skat
    for x in $(seq $No1 $No2) X
    do
        chr="chr"$x
        echo ${sample} " go to " ${chr}
        ${PLINK} --vcf ${oPath}/${sample_vcf}.${chr}.vcf --threads 128 --pheno  ${HOME}/${sample_vcf}.fam  --mpheno 4 --allow-no-sex --make-bed --out ${oPath}/${sample_vcf}.${chr}
        awk 'FNR>1{split($8,b,";");split(b[1],a,"="); if ($1!~"#") print a[2]"\t"$3}' ${oPath}/${sample_vcf}.${chr}.vcf >  ${oPath}/${sample_vcf}.${chr}.SKAT
        Rscript ${HOME}/R_SKAT.R -i ${oPath}/${sample_vcf}.${chr} -c ${HOME}/SKAT_188R.cov -o ${oPath}/skat/${sample_vcf}.${chr}.vcf.R.skata_o 
    done
    
    echo -e "${sample} skat-o endTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
fi

#------------------------------------------------------------------------------------------------------------------------------
## PLINK QTL variant association test
#------------------------------------------------------------------------------------------------------------------------------

if [ -n "${con}" ] && ([ "${con}" == "qtl" ]|| [ "${con}" == "#all" ]); then
    mkdir ${oPath}/qtl_sgl
    echo "############################################# PLINK QTL variant association test ##################################"
    echo -e "${sample} plink-qtl startTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    for n in $(seq $No1 $No2) X
    do
        chr="chr"$n
        echo ${sample} " go to " ${chr}
        
        N1=$(awk 'FNR==1{print NF}' ${HOME}/${sample_vcf}.protein)
        N2=$((N1-2))
        
        ${PLINK} --vcf ${oPath}/${sample_vcf}.${chr}.vcf --keep ${HOME}/${sample_vcf}.protein --make-bed --out ${oPath}/${sample_vcf}.${chr}.vcf.gene.soma.R
        
        for x in $(seq 1 $N2)
        do       
            ${PLINK} --bfile ${oPath}/${sample_vcf}.${chr}.vcf.gene.soma.R --pheno ${HOME}/${sample_vcf}.protein --mpheno ${x}  --covar ${HOME}/${sample_vcf}.cov --covar-number  2 --threads 128  --allow-no-sex --linear  --assoc perm --pfilter 0.05 --out ${oPath}/qtl_new/${sample_vcf}.${chr}.R.${x} 
        done
    done
    
    rm ${oPath}/qtl_sgl/*.log
    rm ${oPath}/qtl_sgl/*sex
    
    echo -e "${sample} plink-qtl endTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
fi

#------------------------------------------------------------------------------------------------------------------------------
## PLINK QTL gene burden test
#------------------------------------------------------------------------------------------------------------------------------

if [ -n "${con}" ] && ([ "${con}" == "qtl_gb" ]|| [ "${con}" == "all" ]); then
    mkdir ${oPath}/qtl_gb
    echo "############################################# PLINK QTL gene burden test ##################################"
    echo -e "${sample} plink-qtl gene burden startTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
    
    for n in $(seq $No1 $No2) X
    do
        chr="chr"$n
        echo ${sample} " go to " ${chr}
        
        No=$(ls ${oPath}/${sample_vcf}.${chr}.vcf.geneset_* | wc -l)
        No3=$((No-1))
        N1=$(awk 'FNR==1{print NF}' ${HOME}/${sample_vcf}.protein)
        N2=$((N1-2))
        
        ${PLINK} --vcf ${oPath}/${sample_vcf}.${chr}.vcf --keep ${HOME}/${sample_vcf}.protein --recode --make-bed --out ${oPath}/${sample_vcf}.${chr}.vcf.gene.qtl.R
        
        for x in $(seq 1 $N2)
        do
            for i in $(seq 0 $No3)
            do
                ${PLINK} --bfile ${oPath}/${sample_vcf}.${chr}.vcf.gene.qtl.R --pheno ${HOME}/${sample_vcf}.protein --mpheno ${x} --threads 128     --covar ${HOME}/${sample_vcf}.cov --covar-number  2 --allow-no-sex --linear --set ${oPath}/${sample_vcf}.${chr}.vcf.geneset_${i}  --assoc perm set-test --set-p 0.05  --pfilter 0.05 --out ${oPath}/qtl_gb/${sample_vcf}.${chr}.vcf.set.R.${x}_${i} 
            done
        done
    done
    
    rm ${oPath}/qtl_gb/*.log
    rm ${oPath}/qtl_gb/*sex
    echo -e "${sample} plink-qtl gene burden endTime::::::::\c" >>${oPath}/plink_${date}.log && date >>${oPath}/plink_${date}.log
fi
