## Captivity and the Co-diversification of Great Ape Microbiomes
### Authors: Alex H. Nishida & Howard Ochman


#### This manuscript is currently is preparation.

#### About this repository
This repository contains scripts for processing raw fastq data as well as scripts for analyzing data to generate figures and tables for the manuscript.

#### Original fastq files:
The original Fastq files were submitted to the NCBI sequence read archive under BioProject accession numbers PRJNA692991, PRJNA693013.

#### Published fastq data:
Published fastq metadata are listed in Table S1

#### File structure
This study analyzes two seperate datasets - 16S and gyrb amplicon datasets.

##### Required programs and dependencies to process 16S amplicon data

###### V4_demultiplex.sh: demultiplex(v1.0.1), repair.sh from bbmap (v38.70), barcode-splitter(v0.18.6)
###### V4_cutadapt.sh: cutadapt(v2.5)
###### DADA2_V4.R, DADA2_single.R, DADA2_paired.R: DADA2(v1.16.0), phyloseq(v1.32.0), tidyverse(v1.3.0)
###### merge_DADA_to_Phyloseq_V4.R
Rpackages: DADA2(v1.16.0), phyloseq(v1.32.0), genefilter(v1.70.0), tidyverse(v1.3.0)
seqinr(v3.6-1), ape(v5.4-1), phytools(v0.7-70)
Programs: FastTree(v2.1.9), mafft(v7.309)

##### Required programs and dependencies to process gyrB amplicon data
The gyrb dataset uses all of the programs listed above plus some additional programs to assign taxonomy

###### create_Bacteroidetes_GTDBTK_ref.ipynb: 
jupyter notebook
python packages: Biopython(v1.77)
programs=prodigal(v2.6.3),hmmer(v3.3)

###### moeller_sup_codiv_clades.ipynb: python packages:Biopython(v1.77),ete3(v3.1.1)
###### gyrb_cutadapt.sh: cutadapt(v2.5)
###### DADA2_gyrb_Bt.R, DADA2_single.R: DADA2(v1.16.0), phyloseq(v1.32.0), tidyverse(v1.3.0)

###### filter_gyrb_seqs_from_metagenomic_samples.ipynb
###### generate_metadata_for_metagenomic_samples.ipynb
###### scripts/processing/gyrb/merge_amplicon_metagenomic_datasets.R
###### idTaxa.R
###### blastp_filter_ASVs.sh 

##### Required programs and dependencies used by scripts in the analyses folder:
-Uses data in the inputs folder to perform all analyses, generatefigures and tables
