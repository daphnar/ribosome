# This package converts CRAM input files to fastq files containing rDNA reads which are entered as input to RiboVAn to extract variant frequencies
## This Applet runs on DNAnexus cloud on the UK Biobank dataset  

![RiboVAn Illustration](https://github.com/daphnar/ribosome/blob/main/rRNA/RiboVAn/RiboVAn_illustration.png)

***
This applet runs on a chunk of 100 files, and should be run in parallel for 500 chunks to get the full UK Biobank cohort  
The main script src/code.sh downloads the human reference genome, prepares the environment, and runs single_run.sh in parallel
***
### Example for running rdna_var_from_cram on the first 100 individuals of the UK Biobank (chunk 0)
> dx run -ichunk=0 -y rdna_var_from_cram --priority low --destination rDNA\ Variations:output_by_mount
