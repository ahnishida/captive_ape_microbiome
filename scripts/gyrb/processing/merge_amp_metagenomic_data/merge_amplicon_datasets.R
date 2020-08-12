library(dada2); packageVersion("dada2")
library(phyloseq)
library(tidyverse)
library(phytools)
library(seqinr)
library("genefilter")
setwd('/Volumes/AHN/captive_ape_microbiome')

#this merges amplicon datasets and filters out reads that aren;t assign. 
#filtered ASVs are then merged with the gtdbtk reference fasta and used to blastn query metagenomic datasets

#load assign taxonomy script
source('scripts/gyrb/processing/merge_amp_metagenomic_data/assign_taxonomy.R')

#define inputs/outputs
outpath='results/gyrb/processing/gyrb_amp_datasets/'
gyrb_ref_fasta = 'results/gyrb/processing/ref_gyrb_gtdbtk/gyrb_fastas/gtdbtk_gyrb_Bt_amplicon.fasta'
dir.create(outpath, recursive=TRUE)

#import data from separate gyrb amplicon runs and human metagenomic data
print('import gyrb_moeller_wild dataset')
gyrb_moeller_wild_path <- "results/gyrb/processing/gyrb_moeller_wild/DADA2" 
gyrb_moeller_wild.seqtab.nochim <- readRDS(file.path(gyrb_moeller_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_moeller_wild.seqtab.nochim )
print('import gyrb_nishida_captive_wild dataset')
gyrb_nishida_captive_wild_path <- "results/gyrb/processing/gyrb_nishida_captive_wild/DADA2" 
gyrb_nishida_captive_wild.seqtab.nochim <- readRDS(file.path(gyrb_nishida_captive_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_nishida_captive_wild.seqtab.nochim)
print('import gyrb_nishida_projectchimps dataset')
gyrb_nishida_projectchimps_path <- "results/gyrb/processing/gyrb_nishida_projectchimps/DADA2" 
gyrb_nishida_projectchimps.seqtab.nochim <- readRDS(file.path(gyrb_nishida_projectchimps_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_nishida_projectchimps.seqtab.nochim )

#merge ASV tables
gyrb.asv <- mergeSequenceTables(gyrb_moeller_wild.seqtab.nochim,
                                gyrb_nishida_captive_wild.seqtab.nochim,
                                gyrb_nishida_projectchimps.seqtab.nochim)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(gyrb.asv)
asv_headers <- vector(dim(gyrb.asv)[2], mode="character")
for (i in 1:dim(gyrb.asv)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(outpath,"ASVs_all.fasta"))

#inputs
asv_fasta_file = file.path(outpath,"ASVs_all.fasta")
asv_faa_file = file.path(outpath,"ASVs_all.faa")
paste0('transeq -frame 2 -sequence ',asv_fasta_file,' -outseq ',asv_faa_file)

outfolder = file.path(outpath,'assign_taxonomy')
pident_cutoff <- 50
length_cutoff <- .8
assign_taxonomy(asv_fasta_file,asv_faa_file,outfolder,pident_cutoff,length_cutoff) 

tax_table <- read.table(paste0(outfolder,'/ASVs_taxonomy.txt'),header = TRUE)
tax_table_ASVs <- tax_table$ASV

#output select ASV fasta
(select_asv_indices = which(asv_headers %in% paste0(">", tax_table_ASVs)))
(select_asv_headers = asv_headers[select_asv_indices])
(select_asv_seqs = asv_seqs[select_asv_indices])
(select_asv_fasta <-  c(rbind(select_asv_headers, select_asv_seqs)))
write(select_asv_fasta, file.path(outpath,"ASVs_filtered.fasta"))

system(paste0('cat ',gyrb_ref_fasta,' ',file.path(outpath,"ASVs_filtered.fasta"),' > ',file.path(outpath,"Bt_ASV_gtdbtk.fasta")))
#cp Bt_ASV_gtdbtk.fasta into metagenomic_samples folder