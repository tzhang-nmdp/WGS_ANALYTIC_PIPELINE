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

* `2. PLINK gene burden test;`

* `3. R conditional logistic regression test;`

* `4. R SKAT-O test;`

* `5. PLINK QTL test;`

* `6. PLINK QTL gene burden test.`

  
#### PART III. ASSOCIATION TEST OUTPUT CONTENT CHECK

* '1. Output for PLINK logstic/Chisq test: ***.assoc.logistic.perm /***.assoc.perm'
<img src="https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_logistic_output.png" width="800">
* '2. Output for PLINK gene burden test: ***.assoc.set.perm'
<img src="https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_logistic_gene_burden_output.png" width="800">
* '3. Output for R conditional logistic regression test: ***.clr'
<img src="https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/CLR__output.png" width="800">
* '4. Output for R SKAT-O test'
<img src="https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/SKAT_O_output.png" width="800">
* '5. Output for PLINK QTL test'
<img src="https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_QTL_output.png" width="800">
* '6. Output for PLINK QTL gene burden test'
<img src="https://github.com/taozhang2019/WGS_ANALYTIC_PIPELINE/blob/master/example/PLINK_QTL_gene_burden_output.png" width="800">
