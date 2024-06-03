#!/usr/bin/env Rscript

#set up workdir

library(survival)
library(optparse)

print("#optional arguments")
option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"),make_option(c("-d", "--head_file"),type="character",help="header file", default=NA,metavar="filename"),make_option(c("-s","--sample_file"),type="character",help="sample_file",default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
input<-as.character(opt$input_file)
head_f<-as.character(opt$head_file)
sample_info<-as.character(opt$sample_file)
output<-as.character(opt$output_file)

cwd=dirname('output')
setwd(cwd)

print("starting...")
Sys.time()

print("# read ID pair group information")
id<-read.table(sample_info,header=T, stringsAsFactors = F)
id_T<-t(id)

print("# read variant / sample information from vcf file")
var_sample<-read.table(input,header=F,sep="\t",stringsAsFactors = F)

print("# read sample header for vcf table ")
sam_h<-read.table(head_f,header=F,stringsAsFactors = F)

colnames(var_sample)<-sam_h
v_len<-length(var_sample$ID)
s_len<-length(var_sample[1,])

print("# transform genotype information")
varid_list<-c()
pvalue_list<-c()
geno_factor_replace<-function(n)
{
if (grepl('0/0',n))
{n<-as.factor('0')}
else if (grepl('./.',n))
{n<-as.factor('0')}
else
{n<-as.factor('1')}
}


for ( i in 1:v_len)
{
print("# integrate genotype phenotype information")
    varid_list[i]<-var_sample[i,3]
    var_tmp0<-t(var_sample[i,c(10:s_len)])
    var_tmp<-as.data.frame(var_tmp0)
    var_tmp[,2]<-row.names(var_tmp)
    colnames(var_tmp)<-c('V1','V2') 
    id_var<-as.data.frame(merge(id,var_tmp,by.x='IID',by.y='V2'))
    ln<-dim(id_var)[2]
    colnames(id_var)[ln]<-'geno'  
    id_var[,ln]<-sapply(id_var[,ln],as.character)         
    id_len<-dim(id_var)[1]  
    for ( n in 1:id_len)
    {
    if (grepl('0/0',id_var[n,ln]))
        {
        id_var[n,ln]<-0}
    else
        {
        id_var[n,ln]<-1}
        }        
    id_var[,ln]<-sapply(id_var[,ln],as.factor)
    sm_len<-length(id_var$geno)
print("# clogit test running")    
    if (sum(as.numeric(as.character(id_var$geno)))==0)  
        {pvalue_list[i]<-1}
    else if (sum(as.numeric(as.character(id_var$geno)))==sm_len)
        {pvalue_list[i]<-1}           
    else
        {
    var_CLR<-clogit(outcome~geno+strata(COV1),method="approximate", data=id_var)
    warnings()
    st_tmp<-summary(var_CLR)
    pvalue_list[i]<-st_tmp$coefficients[5]
        }
}

print("# clogit test output")  
var_list<-cbind(varid_list,pvalue_list)
write.table(var_list,output,quote=FALSE,row.names = FALSE)
print("Ending...")
Sys.time()
