import sys,os,re
import multiprocessing as mp

# basic arguments
vcftyper=open(sys.argv[1],'r')
con=str(sys.argv[2])

# define function for correction of somatic genotype and info field 
def vcf_id_somatic(line):
    item=line.strip().split('\t')
    if item[0].startswith('##'):        
        return "NONE"    
    elif item[0].startswith('#CHROM'):
        line='\t'.join(item)
        return line
    else:    
        if ':' not in item[2]:   
            item[2]=item[0]+':'+item[1]+':'+item[3]+':'+item[4]          
        l_len=len(item)
        for i in range(9,l_len):
            if item[i]=='./.:.:.:.:.' or item[i]=='.|.:.:.:.:.':
                item[i]='0/0:0,0'        
            geno=item[i].split(':')[0]
            geno_n=re.sub('\|','/',geno)
            allele=(item[i].split(':')[1]).split(',')
            for k,v in enumerate(allele):
                if v=='.':
                    allele[k]=0            
            m1=sorted(allele)[-1]
            ref=int(allele[0])
            for n in range(0,len(allele)):
                if allele[n]==m1:
                    alt=int(allele[n])
                    m_1=n               
            if m_1!=0 and alt<5: #or ref>=20 
                item[i]="0/0"+':'+','.join(str(w) for w in allele)   
            elif geno.startswith('./.') and alt>=5: #and ref>=20 
                item[i]="0/"+str(m_1)+':'+str(ref)+','+str(alt)                                                    
            else:
                item[i]=geno_n +':'+','.join(str(w) for w in allele)                                      
        info=item[7].split('|')
        m=(len(info)-1)/15+1
        gene=''
        for i in range(1,m):
            g=int(i*15-12) 
            gene=info[g]
            v_type=info[g-2]            
            if gene!='':
                break  
        if gene=='':
            gene="NNNNN"  
        item[7]='GENE='+gene+';'+'ANN='+v_type
        line='\t'.join(item)
        return line

# run vcf transformation function                 
pool=mp.Pool(processes=32)
if con=="somatic":
    output=[pool.apply(vcf_id_somatic,args=(line,)) for line in vcftyper]
    for line in output:
        if not line.startswith('NONE'):
            print(line)        
