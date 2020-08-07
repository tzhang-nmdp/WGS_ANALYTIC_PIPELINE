#!/usr/bin/env Rscript
setwd("/home/tzhang")

library(SKAT)
library(optparse)

print("#optional arguments")
option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"),make_option(c("-c", "--cov_file"),type="character",help="covariate file", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
input<-as.character(opt$input_file)
#sample_info<-as.character(opt$sample_file)
#outputc<-paste(as.character(opt$output_file),'.c',sep='')
outputo<-paste(as.character(opt$output_file),'.o',sep='')
print("starting...")
Sys.time()

cwd=dirname('output')
setwd(cwd)

print("#prepare the data input in plink format and define the SSD parameter")
File.Fam<- paste(input,'.fam',sep='')
File.Bim<- paste(input,'.bim',sep='')
File.Bed<- paste(input,'.bed',sep='')
File.Set<- paste(input,'.SKAT',sep='')
File.SSD<- paste(input,'.SSD',sep='')
File.Info<-paste(input,'.SSD.Info',sep='')

print("#initiate the SSD setting")
FAM<-Read_Plink_FAM(File.Fam,Is.binary = T)
FAM[is.na(FAM)]<-1
print(FAM)
#FAM[FAM[,6]==0,6]<-1
Generate_SSD_SetID(File.Bed,File.Bim,File.Fam,File.Set,File.SSD,File.Info)
print("ok")
SSD.INFO<-Open_SSD(File.SSD,File.Info)
#y<-FAM$Phenotype

print("#start covariate")
File.Cov<-as.character(opt$cov_file)
FAM_Cov<-Read_Plink_FAM_Cov(File.Fam, File.Cov, cov_header=TRUE, Is.binary=TRUE)
X1 = FAM_Cov$outcome
# supress the covariates due to heavy computation burden
X2 = FAM_Cov$COV1
#X3 = FAM_Cov$COV3
#X4 = FAM_Cov$COV4
#X5 = FAM_Cov$COV5
#X6 = FAM_Cov$COV6
#X7 = FAM_Cov$COV7
#X8 = FAM_Cov$COV8
#X9 = FAM_Cov$COV9
#X10 = FAM_Cov$COV10
#X11 = FAM_Cov$COV11
#X12 = FAM_Cov$COV12
#X13 = FAM_Cov$COV13
#y<-FAM_Cov$PHE

print(input)
print("#start the KAT-O test with resampling adjustment")
#obj<-SKAT_Null_Model(y ~ X1, out_type="C", n.Resampling=1000, type.Resampling="bootstrap")
#outc.skatc<-SKAT_CommonRare(SSD.INFO, obj, r.corr.rare=1, r.corr.common=1)
obj<-SKAT_Null_Model(Phenotype ~ X1+X2, out_type = "D",data=FAM_Cov) #+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12
outo.skato<-SKATBinary.SSD.All(SSD.INFO,obj,method="SKATO") #, N.Resampling=2 *10^6, seednum=100, epsilon=10^-6, SetID=NULL)

print("#SKAT-O output saving")
warnings() 
write.table(outo.skato$result,outputo,quote=FALSE)
print("Ending...")
Sys.time()
