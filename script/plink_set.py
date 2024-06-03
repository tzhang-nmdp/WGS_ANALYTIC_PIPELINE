# prepare gene_set file for PLINK gene burden test and SKAT-O test
import sys

# basic parameters
input_vcf=open(sys.argv[1],'r')
ab_vcf=open(str(sys.argv[2])+'.ab','w')
var_list=[]
gene_dict={}
info_dict={}
i=n=0
try:
    filt_vcf=open(sys.argv[3],'r')
    for line in filt_vcf:
        item=line.strip().split('\t')  
        var_list.append(item[0])
except:
    print("no filtering")
    
# loop all variant gene information and collapse into gene sets
for line in input_vcf:
    item=line.strip().split('\t')  
    if item[0].startswith('#') or item[0].startswith('CHROM') or len(item)<1:
        ab_vcf.write(line)    
        continue      
    var_id=item[2]
    if var_id in var_list:
        ab_vcf.write(line)
        continue
    info=item[7].split(';')
    for inf in info:
        key=inf.split('=')[0]
        value=inf.split('=')[1]            
        info_dict[key]=value
    gene=info_dict['GENE']
    if gene=='NNNNN':
        ab_vcf.write(line)    
        continue
    if gene_dict.get(gene):
        gene_dict[gene].append(var_id)
    else:
        gene_dict.update({gene:[var_id]})
  
# print out gene sets with size <280 (memory limit)        
for key, value in gene_dict.items():
    if (i % 280 ==0):
        out_vcf_name=str(sys.argv[2])+'_'+str(n)
        out_vcf=open(out_vcf_name,'w')
        n=n+1
    out_vcf.write(key+'\n'+'\n'.join(value)+'\nEND\n\n')        
    i=i+1
    
