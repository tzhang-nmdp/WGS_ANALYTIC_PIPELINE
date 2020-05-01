# for sample_id fix
import sys
input_vcf=open(sys.argv[1],'r')
sample_file=open(sys.argv[2],'r')
output=open(str(sys.argv[1])+'.cr','w')
sample_list=[]

# generate sample list
for line in sample_file:
    item=line.strip().split('\t')  
    sample_list.append(item[0])  
        
# fix sample id in vcf        
for line in input_vcf:
    item=line.strip().split('\t')  
    if item[0].startswith('#CHROM'):
        out_vcf.write('\t'.join(item[0:9])+'\t'+'\t'.join(sample_list)+'\n')         
    else:
        out_vcf.write('\t'.join(item)+'\n')    
        continue
