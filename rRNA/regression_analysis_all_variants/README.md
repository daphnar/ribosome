# This package runs regression models, either linear or logistic regressions, between all heritable rRNA variants and selected phenotypes
## This Applet runs on DNAnexus cloud on the UK Biobank dataset  
This applet takes as input two parameters, the type of regression model (linear or logistic), and a chunk ID which corresponds to a subset of phenotypes

### Example usage for regression_type logistic and chunk 0
> dx run -ichunk=0 -iregression_type=logistic -y regression_analysis --priority low --destination rDNA\ Variations:regression_analysis --instance-type mem1_hdd1_v2_x16


![GWAS example](https://github.com/daphnar/rRNA/blob/main/regression_analysis_all_variants/regression.png)


