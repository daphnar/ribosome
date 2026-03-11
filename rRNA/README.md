# This repository includes the computational tools published in:  
##  [Rothschild et al. Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100629)
## [Rothschild et al. MedRxiv 2025](https://doi.org/10.1101/2025.09.02.25334953)

In [Rothschild et al. Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100629) I released the Reference Gap Alignment (RGA) algorithm suitible for sequence variant calling in paralog genes.  
It's implemetation and detailed instructions are found under the RGA folder.  

In  [Rothschild et al. MedRxiv 2025](https://doi.org/10.1101/2025.09.02.25334953) I released the Ribosome Variation Analysis (RiboVAn) pipeline which offers mapping short-reads to the nucleotide resolution rRNA atlas ([Rothschild et al. Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100629)).  
It's implemetation and detailed instructions are found under the RiboVAn folder.  

Additionally, you will find here 3 applets that are suitible for running over the DNAnexus cloud platform on the UK Biobank cohort.  
These applets, found as three folders in this repository, were used in [Rothschild et al. MedRxiv 2025](https://doi.org/10.1101/2025.09.02.25334953) for studying the affect of rRNA variants on human health in the UK Biobank dataset.  

Here is an illustration and motivation for studying the rRNA genes and their sequence variations.  

In [Rothschild et al. Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100629) we created an atlas of variants (panel A in the illustration figure) using the RGA method (implemented here) and RIBO-RT which is an experimental method for sequencing the 18S and 28S rRNA from translating ribosomes (described [here](https://doi.org/10.1016/j.xgen.2024.100629)).  

With this atlas (panel A) we asked if rRNA variants affect human health (panel B) using the whole genome sequencing and rich phenotypic measurements available in the UK Biobank dataset.
For this I wrote the RiboVAn pipeline (panel C), which outputs single nucleotide rRNA variants frequencies which I associated with phenotypes (panel C).  

![rRNA_variant_illustration](https://github.com/daphnar/ribosome/rRNA/rRNA_variant_Illustration.png)

Please cite [Rothschild et al. Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100629) if you are using RGA.  
Please cite [Rothschild et al. MedRxiv 2025](https://doi.org/10.1101/2025.09.02.25334953) if you are using RiboVAn.

For technical support write daphna@stanford.edu
