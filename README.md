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
Scripts in the processing folder:\n
-Perform QC filtering, remove primers, error correction on raw fastq\n
-Generate ASV table, rarefy, assign taxomony, generate phylogeny,\n correct metadata into the inputs fold\n
Scripts in the analysis folder:\n
-Uses data in the inputs folder to perform all analyses, generate\n figures and tables
