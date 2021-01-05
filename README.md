## Captivity and the Co-diversification of Great Ape Microbiomes
### Authors: Alex H. Nishida & Howard Ochman


#### This manuscript is currently is preparation. 

#### About this repository
This repository contains scripts for processing raw fastq data as well as scripts for analyzing data to generate figures and tables for the manuscript.

#### Original fastq files:
The original Fastq files were submitted to the NCBI sequence read archive under BioProject accession number [Insert accession number here].

#### Published fastq data: 
Published fastq metadata are listed in Table 1

#### File structure
This study analyzes two seperate datasets - 16S and gyrb amplicon datasets. 
Scripts in the processing folder:
-Perform QC filtering, remove primers, error correction on raw fastq
-Generate ASV table, rarefy, assign taxomony, generate phylogeny, correct metadata into the inputs fold
Scripts in the analysis folder:
-Uses data in the inputs folder to perform all analyses, generate figures and tables

