library(dada2); packageVersion("dada2")
library(phyloseq)
library(DECIPHER)
library(tidyverse)
library(phytools)
library(seqinr)
library(Biostrings)
#library("genefilter")
setwd('/Volumes/AHN/captive_ape_microbiome')

#COMPARE MAFFT VS TRANSALIGN

setwd('results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2')

mafft_aln <- Biostrings::readDNAStringSet(file.path("GemRef_AlignMafft","ASVs_filtered_ref_aln_trim.fasta"))
aligntrans_aln <- Biostrings::readDNAStringSet(file.path("GemRef_AlignTrans","ASVs_filtered_ref_aln_trim.fasta"))

mafft_aln_no_gaps <- as.character(RemoveGaps(mafft_aln,removeGaps = "all"))
aligntrans_aln_no_gaps <- as.character(RemoveGaps(aligntrans_aln,removeGaps = "all"))

common_seqs <- mafft_aln_no_gaps[mafft_aln_no_gaps %in% aligntrans_aln_no_gaps]
diff_mafft <- mafft_aln_no_gaps[!mafft_aln_no_gaps %in% aligntrans_aln_no_gaps]
diff_mafft['ASV_18']
diff_transalign <- aligntrans_aln_no_gaps[!aligntrans_aln_no_gaps %in% mafft_aln_no_gaps]
diff_transalign['ASV_18']

fasta <- Biostrings::readDNAStringSet(file.path("GemRef_AlignMafft","ASVs_filtered.fasta"))
fasta <- DNAStringSet(fasta,start=2)
names(fasta[!fasta %in% mafft_aln_no_gaps])
names(fasta[!fasta %in% aligntrans_aln_no_gaps])

#seqs dont match ASVs bc ASVs were trimmed to the region that ASV_1 aligns to
#check out the ASVs_filtered_ref_aln.fasta file
