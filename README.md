## Captivity and the Co-diversification of Great Ape Microbiomes
### Authors: Alex H. Nishida & Howard Ochman


#### This manuscript is currently is under-consideration.

#### About this repository
This repository contains scripts for processing raw fastq data as well as scripts for analyzing data to generate figures and tables for the manuscript.

#### Original fastq files:
The original Fastq files were submitted to the NCBI sequence read archive under BioProject accession numbers PRJNA692991, PRJNA693013.

#### Published fastq data:
Published fastq metadata are listed in Table S1

#### File structure
###### This study analyzes two seperate datasets - 16S and gyrb amplicon datasets.

###### Scripts are divided into two folders processing and analyses
###### Processing scripts convert raw fastq to ASV table with associated taxonomy, metadata and phylogeny
###### Analyses scripts use ASV table with associated taxonomy, metadata and phylogeny in the inputs folder to generate all figures and run all statistical tests
###### 
###### Dependencies can be found in the dependencies.md file

###### To reproduce all figures and statistics reported in the manuscript
1) Check the dependencies.md file to obtain python modules and R libraries necessary for running analyses scripts
2) git clone https://github.com/ahnishida/captive_ape_microbiome.git
3) gunzip results/gyrb/inputs/physeq_Bacteroidales_asv_tab.txt.gz
4) Edit working directory in these scripts: 16S_analyses_all.R /path/to/captive_ape_microbiome
5)  