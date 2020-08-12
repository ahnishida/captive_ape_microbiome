library(dada2); packageVersion("dada2")
library("genefilter")
library(tidyverse)
library(phyloseq)
library(seqinr)
library(ape)
library(phytools)

setwd('/Volumes/AHN/captive_ape_microbiome')
inpath = 'results/16s/processing/16s_all/'
outpath = 'results/16s/inputs/'
physeq <- readRDS(file.path(inpath,'phyloseq_rare10000.rds'))
tax_table_tab <- as.data.frame(tax_table(physeq)) %>% 
  rownames_to_column(var='ASV')
write.table(tax_table_tab,file.path(outpath,'ASVs_taxonomy.txt'),quote=F,row.names=F,sep='\t')

metadata_tab <- as.data.frame(sample_data(physeq)) 
metadata_tab$X.SampleID <- rownames(metadata_tab)
write.table(metadata_tab,file.path(outpath,'16S_metadata.txt'),quote=F,row.names=F,sep='\t')

otu_table <- as.data.frame(otu_table(physeq)) 
write.table(otu_table,file.path(outpath,'ASV_tab.txt'),quote=F,row.names=T,col.names=NA,sep='\t')

