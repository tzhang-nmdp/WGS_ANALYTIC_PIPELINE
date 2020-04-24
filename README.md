# WGS_ANALYTIC_PIPELINE

It's WGS analytic pipeline with two parts:

<div class="text-blue">
PART I. A series of VCF data transformations: <a href="#" class="text-inherit">including the link</a>
</div>
PART I. A series of VCF data transformations:
1. snpEff VCF annotation;
2. bcftools VCF merge;
3. VCF data transformation on genotype/annotation fields.

PART II. A set of variant level or gene level association tests:
1. VCF data transformation;
2. PLINK logstic/Chisq test;
3. R conditional logistic regression test;
4. PLINK gene burden test;
5. R SKAT-O test;
6. PLINK QTL test;
7. PLINK QTL gene burden test.
