library(dada2); packageVersion("dada2");
library(Rcpp)
library(phyloseq)
library(tidyverse)
setwd('/Volumes/AHN/captive_ape_microbiome')
source('scripts/DADA2_QC/DADA2_paired.R')
source('scripts/DADA2_QC/DADA2_single.R')

metadata_file='metadata/metadata_all_samples_16S.txt'
metadata <- read.csv(metadata_file,sep='\t')
getid <- function(old_id){
  new_id= metadata$X.SampleID[metadata$Old_SampleID==old_id]
  return(new_id)
}

###CLAYTON###
#define input fastqs, sample names, outdir
rawdata_path <- "data/16S_clayton_captive/raw_fastq" 
outdir <- "results/16S_clayton_captive/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern="_R1_001.fastq.gz", full.names = TRUE)))
(reverse_reads <- sort(list.files(rawdata_path, pattern="_R2_001.fastq.gz", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1))
(sample.names <- as.character(unlist(lapply(sample.names,getid)))) #get new sample names 
#View read quality
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
trunc_parameters = c(145,145) #define truncation parameters
#run DADA2 pipeline
DADA2_paired(forward_reads,reverse_reads, sample.names, outdir, trunc_parameters) 

###GOODRICH####
#define input fastqs, sample names, outdir
rawdata_path <- "data/16s_goodrich/raw_fastq" 
outdir <- "results/16s_goodrich/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern="_1.fastq.gz", full.names = TRUE)))
(reverse_reads <- sort(list.files(rawdata_path, pattern="_2.fastq.gz", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1))
(sample.names <- as.character(unlist(lapply(sample.names,getid)))) #get new sample names 
#View read quality
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
trunc_parameters = c(145,145) #define truncation parameters
#run DADA2 pipeline
DADA2_paired(forward_reads,reverse_reads, sample.names, outdir, trunc_parameters) 

###NISHIDA_CAPTIVE####
rawdata_path <- "results/16s_nishida_captive/demultiplexed_fastq" 
outdir <- "results/16s_nishida_captive/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern="-read-1.fastq", full.names = TRUE)))
(reverse_reads <- sort(list.files(rawdata_path, pattern="-read-2.fastq", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "-"), `[`, 1))
#View read quality
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
trunc_parameters = c(145,145) #define truncation parameters
#run DADA2 pipeline
DADA2_paired(forward_reads,reverse_reads, sample.names, outdir, trunc_parameters) 

###NISHIDA_PROJECTCHIMPS####
#define input fastqs, sample names, outdir
rawdata_path <- "results/16s_nishida_projectchimps/cutadapt_fastq" 
outdir <-  "results/16s_nishida_projectchimps/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern="_R1_001.fastq", full.names = TRUE)))
(reverse_reads <- sort(list.files(rawdata_path, pattern="_R2_001.fastq", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1))
(sample.names <- as.character(unlist(lapply(sample.names,getid))))
#View read quality
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
trunc_parameters = c(145,145) #define truncation parameters
#run DADA2 pipeline
DADA2_paired(forward_reads,reverse_reads, sample.names, outdir, trunc_parameters) 

###RAYMANN###
#define input fastqs, sample names, outdir
rawdata_path <- "results/16s_raymann_captive/cutadapt_fastq" 
outdir <-  "results/16s_raymann_captive/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern="_R1_001.fastq", full.names = TRUE)))
(reverse_reads <- sort(list.files(rawdata_path, pattern="_R2_001.fastq", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "_R1_001.fastq"), `[`, 1))
(sample.names <- as.character(unlist(lapply(sample.names,getid))))
#View read quality
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
trunc_parameters = c(200, 100) #define truncation parameters
#run DADA2 pipeline
DADA2_paired(forward_reads,reverse_reads, sample.names, outdir, trunc_parameters) 

###SMITS####
#define input fastqs, sample names, outdir
rawdata_path <- "data/16s_smits_2017/raw_fastq"  
outdir <-  "results/16s_smits_2017/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern="_1.fastq.gz", full.names = TRUE)))
(reverse_reads <- sort(list.files(rawdata_path, pattern="_2.fastq.gz", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1))
(sample.names <- as.character(unlist(lapply(sample.names,getid))))
#View read quality
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
trunc_parameters = c(145, 145) #define truncation parameters
#run DADA2 pipeline
DADA2_paired(forward_reads,reverse_reads, sample.names, outdir, trunc_parameters) 

##VANGAY###
rawdata_path <- "results/16s_vangay_2018/cutadapt_fastq"   
outdir <-  "results/16s_vangay_2018/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern="R1.fastq.gz", full.names = TRUE)))
(reverse_reads <- sort(list.files(rawdata_path, pattern="R2.fastq.gz", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), ".R1.fastq.gz"), `[`, 1))
(forward_reads <- forward_reads[sample.names %in% metadata$Old_SampleID]) #accidently downloaded more fastq than intended
(reverse_reads <- reverse_reads[sample.names %in% metadata$Old_SampleID])
(sample.names <- sapply(strsplit(basename(forward_reads), ".R1.fastq.gz"), `[`, 1))
(sample.names <- as.character(unlist(lapply(sample.names,getid))))

#View read quality
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
trunc_parameters = c(145, 145) #define truncation parameters
#run DADA2 pipeline
DADA2_paired(forward_reads,reverse_reads, sample.names, outdir, trunc_parameters) 

##MOELLER_CHIMP###
#define input fastqs, sample names, outdir
rawdata_path <- "results/16s_moeller_wild/chimp/cutadapt_fastq"   
outdir <-  "results/16s_moeller_wild/chimp/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern=".fastq.gz", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "-"), `[`, 1))

#View read quality
plotQualityProfile(forward_reads[1:2])
trunc_parameters = c(250)
#no truncation length specified
#run DADA2 pipeline
DADA2_single(forward_reads, sample.names, outdir, trunc_parameters)

###MOELLER_GORBON###
#define input fastqs, sample names, outdir
rawdata_path <-  "results/16s_moeller_wild/gor_bon/demultiplexed_fastq" 
outdir <-  "results/16s_moeller_wild/gor_bon/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern=".fastq.gz", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), ".fastq.gz"), `[`, 1))

#View read quality
plotQualityProfile(forward_reads[1:2])
trunc_parameters = c(250)
#no truncation length specified
#run DADA2 pipeline
DADA2_single(forward_reads, sample.names, outdir, trunc_parameters)

################


