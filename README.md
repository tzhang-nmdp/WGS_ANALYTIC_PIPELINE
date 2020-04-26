# WGS_ANALYTIC_PIPELINE 

|Pre- requisites|
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

