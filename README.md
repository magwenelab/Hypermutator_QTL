# Hypermutator_QTL

Python functions and jupyter notebooks used in QTL analysis of *Cryptococcus neoformans* hypermutator phenotype in **[Priest et al. 2021](https://www.biorxiv.org/content/10.1101/2021.08.11.455996v1)** 

## Dependencies

**[Python (anaconda) v 3.7.3](https://www.anaconda.com/)**
- Used for analysis and visualization

**[BWA v 0.7.12-r1039](http://bio-bwa.sourceforge.net/)**
- Used to align FASTQ file to an XL280 reference genome

**[Samtools v 1.9](http://www.htslib.org/)**
- Used to generate and filter SAM and BAM files

**[Freebayes v 1.2.0](https://github.com/freebayes/freebayes) haplotype caller**
- Used to detect genetic variants segregating in the mapping population

**[Bamaddrg](https://github.com/ekg/bamaddrg)**
- Used to add read group information to BAM

## Description of notebooks and order of analysis pipeline:

### FASTQ_alignment_and_variant_calling
- A notebook that constructs bash calls for aligning paired-end reads with BWA, SAM to BAM conversion, adding readgroups, and calling freebayes for variant detection

### hypermutatorqtl.py
- A library of dataframes, lists, variables, and functions used across notebooks.

### Phenotype_preprocess
- Transforms and addjusts labels and scles of hypermutator phenotype data.

### Figure_S14_Filter_VCF
- Filters the variant call files across chromosomes and generates supplementary Figure S14 in Priest et al.

### Genotype_postprocess
- Reformats genotpe data and sample names for use in QTL mapping.

### Supplementary_Figure_S5
- Plots the haplotypes across chromosome 3 and 11, generating supplementary Figure S5 in Priest et al.

### QTL_analysis
- Codnducts QTL analysis of hyptermutator phenotype and genotypes in Bt65 x H99 F1 progeny, generating upper panel of Figure 2, middle panle of Figure 2, and Supplementary Figure S4 as seen in Priest et al.

### Supplementary_Figure_S13
- For each sample sequenced from the Bt65 x H99 cross, constructs diagnostic plots of allele, read-depth, allelic read-depth across chromosomes. The plot for progney P25 was used as Supplementary Figure S13 in Priest et al.

### LongCOA_visulization
- Plots the gene model, variants, and amino acids for *CNAG_01836* similar to upper panel of Supplementary Figure S7 in Priest et al. 