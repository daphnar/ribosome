# This package implements the Ribosome Variation Analysis (RiboVAn) pipeline

In [Rothschild et al. Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100629) we created an atlas of variants (Panel A in the illustration figure) that allow mapping short reads to the 18S and 28S with variations.
Here I wrote a pipeline Ribosome Variation Analysis (RiboVAn) for getting single nucleotide variants from short reads that map to the atlas (Panel C) which offers for example associating variants with traits (Panel C).

![RiboVAn Illustration](https://github.com/daphnar/rRNA/blob/main/RiboVAn/RiboVAn_illustration.png)

***
The RiboVAn.py script expect as input parameters mapped sam files to the atlas and returns nucleotide atlas variant frequencies
***
As prerequisits to running RiboVAn.py, you need to map short-reads to the atlas of ES/non-ES resolution - found available to download in [Rothschild et al. Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100629).    
You can also use the ES/non-ES resolution atlas already found as a bowtie2 reference in this repository under [bowtie2_atlas_expand150_ES](../rdna_var_from_cram/resources/home/dnanexus/bowtie2_atlas_expand150_ES)
***
The way that RiboVAn works is by extracting the positions of where variants are in the 18S and 28S using a lookup table - all_es.atlas.position.variants.csv
Then, using the nucleotide atlas files, named atlas_18s.XXX and atlas_28s.XXX variant positions and IDs are replaced with matching sequences at found positions.

## Usage example 
### In this example input files were paired-end short reads and we map both strand reads to the atlas. 
fastq_file1=file_1.fastq.gz
fastq_file2=file_2.fastq.gz
#
bowtie2 -x bowtie2_index -U "$fastq_file1" --score-min 'C,0,-1' \
-S "file_1.ribo.atlas_mapped.sam" --un "file_1.unmapped.sam" \
2> "file_1.mapping_stats.txt"
#
bowtie2 -x bowtie2_index -U "$fastq_file1" --score-min 'C,0,-1' \
-S "file_2.ribo.atlas_mapped.sam" --un "file_1.unmapped.sam" \
2> "file_2.mapping_stats.txt"
#
Then file_2.mapping_stats.txt and file_2.mapping_stats.txt are entered as input files to RiboVAn.py as followed:

python RiboVAn.py file_1.ribo.atlas_mapped.sam file_2.ribo.atlas_mapped.sam output
#
The output result will contain nucleodite atlas variant frequencies.
