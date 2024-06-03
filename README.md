# WGS_ANALYTIC_PIPELINE 
![WGS_WORKFLOW](https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/WGS_WORKFLOW.png)

|Pre- requisites|Package information|
|---------------|---------------|
|[miniconda software]: | https://docs.conda.io/en/latest/miniconda.html#linux-installers |
| [snpEff software]: | http://snpeff.sourceforge.net/SnpEff_manual.html  |
| [bcftools software]: | https://samtools.github.io/bcftools/howtos/index.html  | 
| [PLINK software]: | https://www.cog-genomics.org/plink/1.9/output |
| python |
| R (SURVIVAL, SKAT-O) | 

  
#### PART I. A series of VCF data transformations:

* `1. snpEff VCF annotation;`

* `2. bcftools VCF merge;`

* `3. VCF data transformation on genotype/annotation fields;`

* `4.  VCF data transformation on variant gene set.`


#### PART II. A set of variant level or gene level association tests:

* `1. PLINK logstic/Chisq test;`

* `2. R conditional logistic regression test;`

* `3. PLINK gene burden test;`

* `4. R SKAT-O test;`

* `5. PLINK QTL test;`

* `6. PLINK QTL gene burden test.`

  
#### PART III. ASSOCIATION TEST OUTPUT CONTENT CHECK

* '1. Output for PLINK logstic/Chisq test: ***.assoc.logistic.perm /***.assoc.perm'
<p align="center">
  <img src="[your_relative_path_here](https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_logistic_output.png)" width="350" title="">
</p>
* '2. Output for R conditional logistic regression test: ***.clr'
![CLR_output](https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/CLR_output.png)
* '3. Output for PLINK gene burden test: ***.assoc.set.perm'
![PLINK_logistic_gene_burden_output](https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_logistic_gene_burden_output.png) 
* '4. Output for R SKAT-O test'
![SKAT-O_output](https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/SKAT-O_output.png)
* '5. Output for PLINK QTL test'
![PLINK_QTL_output](https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_QTL_output.png) 
* '6. Output for PLINK QTL gene burden test'
![PLINK_QTL_gene_burden_output](https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_QTL_gene_burden_output.png)

