library('dada2'); packageVersion("dada2");
library(Rcpp)
library(phyloseq)
library(tidyverse)
setwd('/Volumes/AHN/captive_ape_microbiome')
#wrapper script to run DADA2_single on 3 gyrB amplicon datasets
#renames samples, selects only Bacteroidaceae samples
#inputs: fastqs after primer trimming
source('scripts/processing/gyrb/amplicon_datasets/DADA2_single.R')
#outputs: asv table where rows are sample names, columns are seqs, and values are read counts

metadata_file='metadata/metadata_gyrb_amplicon_Bt.txt'
metadata <- read.csv(metadata_file,sep='\t')
getid <- function(old_id){
  new_id= metadata$X.SampleID[metadata$Old_SampleID==old_id]
  return(new_id)
}

### MOELLER WILD ###
#define input fastqs, sample names, outdir
rawdata_path <- "results/gyrb/processing/gyrb_moeller_wild/cutadapt_fastq" 
outdir <- "results/gyrb/processing/gyrb_moeller_wild/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern=".fastq", full.names = TRUE)))
(sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1))

#filter out fastqs not present in metadata 
moeller_samples_metadata <- metadata %>% filter(dataset=='gyrb_moeller_wild') 
sample.name_in_metadata <- intersect(sample.names,moeller_samples_metadata$Old_SampleID)
forward_reads_in_metadata <- forward_reads[sample.names %in% sample.name_in_metadata]
sample.names_in_metadata <- as.character(unlist(lapply(sample.names,getid))) #get new sample names 

#now select only Bt samples
(Bt_forward_reads <- forward_reads_in_metadata[stringr::str_detect(sample.names_in_metadata,'Bt')])
(Bt_sample.names <- sample.names_in_metadata[stringr::str_detect(sample.names_in_metadata,'Bt')])
#there is one sample that had the same name as another sample in moeller data but diff barcodes
Bt_sample.names[stringr::str_detect(Bt_sample.names,'dup')] #I'm leaving it in
#View read quality
plotQualityProfile(Bt_forward_reads[1:2])

trunc_parameters = c(250) #define truncation parameters
#run DADA2 pipeline
DADA2_single(Bt_forward_reads, Bt_sample.names, outdir, trunc_parameters) 

###NISHIDA_CAPTIVE####
rawdata_path <- "results/gyrb/processing/gyrb_nishida_captive_wild/cutadapt_fastq" 
outdir <- "results/gyrb/processing/gyrb_nishida_captive_wild/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern=".fastq.gz", full.names = TRUE)))

#now select only Bt samples
(Bt_forward_reads <- forward_reads[stringr::str_detect(forward_reads,'Bt')])
(Bt_sample.names <- sapply(strsplit(basename(Bt_forward_reads), "_"), `[`, 1))
(Bt_sample.names <- as.character(unlist(lapply(Bt_sample.names,getid)))) #get new sample names 

#some samples represent an aggregate of cultured strains from a host species and location
Bt_sample.names[stringr::str_detect(sample.names,'strain')] 
#View read quality
plotQualityProfile(Bt_forward_reads[1:10])

trunc_parameters = c(250) #define truncation parameters
#run DADA2 pipeline
DADA2_single(Bt_forward_reads, Bt_sample.names, outdir, trunc_parameters) 

###NISHIDA_PROJECTCHIMPS####
#define input fastqs, sample names, outdir
rawdata_path <- "results/gyrb/processing/gyrb_nishida_projectchimps/cutadapt_fastq" 
outdir <-  "results/gyrb/processing/gyrb_nishida_projectchimps/DADA2" 
(forward_reads <- sort(list.files(rawdata_path, pattern=".fastq.gz", full.names = TRUE)))
#now select only Bt samples
(Bt_forward_reads <- forward_reads[stringr::str_detect(forward_reads,'Bt')])
(Bt_sample.names <- sapply(strsplit(basename(Bt_forward_reads), "_"), `[`, 1))
(Bt_sample.names <- as.character(unlist(lapply(Bt_sample.names,getid)))) #get new sample names 

#View read quality
plotQualityProfile(Bt_forward_reads[1:2])

trunc_parameters = c(250) #define truncation parameters
#run DADA2 pipeline
DADA2_single(Bt_forward_reads, Bt_sample.names, outdir, trunc_parameters) 

