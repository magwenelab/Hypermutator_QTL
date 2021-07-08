# Hypermutator_QTL

Python functions and notebooks used in QTL analysis of C. neoformans hypermutator phenotype in Priest et al. 2021

## Dependencies
BWA
SAMTOOLS
BAMADDRG
FREEBAYES
PYTHON

## Description of notebooks and order of analysis pipeline:

FASTQ_alignment_and_variant_calling.ipynb
- A notebook that constructs calls for aligning paired-end reads with BWA, SAM to BAM conversion, adding readgroups, and calling freebayes for variant detection

hypermutatorqtl.py
- A library of dataframes, lists, variables, and functions used across notebooks.

Phenotype_preprocess.ipynb
- Transforms and addjusts labels and scles of hypermutator phenotype data.

Figure_S14_Filter_VCF.ipynb
- Filters the variant call files across chromosomes also generates supplementary Figure S14 in Priest et al. 2021.

Genotype_postprocess.ipynb
- Reformats genotpe data and sample names for use in QTL mapping.

Figure_S9.ipynb
- Plots the haplotypes across chromosome 3 and 11, generating supplementary Figure S9 in Priest et al. 2021.

Figure_4_QTL_map.ipynb
- Codnducts QTL analysis of hyptermutator phenotype and genotypes in Bt65 x H99 F1 progeny, generating Figure 4 of Priest et al. 2021

Figure_S13_Genome_diagnostic_plot
- For each sample sequenced from the Bt65 x H99 cross, constructs diagnostic plots of allele, read-depth, allelic read-depth across chromosomes. The plot for progney P25 was used as supplementary Figure S13 in Priest et al. 2021.