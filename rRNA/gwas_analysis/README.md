# This package runs a Genome Wide Association Study (GWAS) between a single rRNA variant and other genomic SNPs
## This applet runs on DNAnexus cloud on the UK Biobank dataset

### It works by associating variant frequencies of the selected variant to genomic SNPs
This is similar to a standard GWAS where the rRNA variant frequencies across the UK Biobank is treated as a phenotype in the GWAS analysis

### Example for the first variant (chunk 0 which is 28s_59_A_SNV) and chromosome 11 
> dx run -ichunk=0 -ichromosome=11 -y gwas_analysis --priority low --destination rDNA\ Variations:GWAS --instance-type mem1_hdd1_v2_x16

![GWAS example](https://github.com/daphnar/ribosome/blob/main/rRNA/gwas_analysis/GWAS.png)
