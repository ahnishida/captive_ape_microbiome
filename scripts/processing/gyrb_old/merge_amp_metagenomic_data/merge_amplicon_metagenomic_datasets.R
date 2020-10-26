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
outpath='results/gyrb/processing/gyrb_amp_meta/'
gyrb_ref_fasta = 'results/gyrb/processing/ref_gyrb_gtdbtk/gyrb_fastas/gtdbtk_gyrb_Bt_amplicon.fasta'
gyrb_ref_full_fasta = 'results/gyrb/processing/ref_gyrb_gtdbtk/gyrb_fastas/gtdbtk_gyrb_Bt.fasta'

dir.create(outpath, recursive=TRUE)

#import data from separate gyrb amplicon runs and human metagenomic data
print('import gyrb_moeller_wild dataset')
gyrb_moeller_wild_path <- "results/gyrb/processing/gyrb_moeller_wild/DADA2" 
gyrb_moeller_wild.seqtab.nochim <- readRDS(file.path(gyrb_moeller_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_moeller_wild.seqtab.nochim)
gyrb_moeller_wild.seqtab.nochim[gyrb_moeller_wild.seqtab.nochim>0] <- 1
print('import gyrb_nishida_captive_wild dataset')
gyrb_nishida_captive_wild_path <- "results/gyrb/processing/gyrb_nishida_captive_wild/DADA2" 
gyrb_nishida_captive_wild.seqtab.nochim <- readRDS(file.path(gyrb_nishida_captive_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_nishida_captive_wild.seqtab.nochim)
gyrb_nishida_captive_wild.seqtab.nochim[gyrb_nishida_captive_wild.seqtab.nochim>0] <- 1
print('import gyrb_nishida_projectchimps dataset')
gyrb_nishida_projectchimps_path <- "results/gyrb/processing/gyrb_nishida_projectchimps/DADA2" 
gyrb_nishida_projectchimps.seqtab.nochim <- readRDS(file.path(gyrb_nishida_projectchimps_path,"seqtab.nochim.rds")) #Opens ASV table
gyrb_nishida_projectchimps.seqtab.nochim[gyrb_nishida_projectchimps.seqtab.nochim>0] <- 1
rownames(gyrb_nishida_projectchimps.seqtab.nochim )
print("import human_metagenomic_data")
metagenomic_path <- "results/gyrb/processing/metagenomic_samples/final_outputs/gyrb_metagenomic_seqtab.txt"
metagenomic_human_reads <- as.matrix(read.table(metagenomic_path,header=T,
                                                row.names=1, check.names=F, sep="\t")) #Opens ASV table
rownames(metagenomic_human_reads)
colnames(metagenomic_human_reads)
rowSums(metagenomic_human_reads)

#merge ASV tables
gyrb.asv <- mergeSequenceTables(gyrb_moeller_wild.seqtab.nochim,
                                gyrb_nishida_captive_wild.seqtab.nochim,
                                gyrb_nishida_projectchimps.seqtab.nochim,
                                metagenomic_human_reads)
dim(gyrb.asv)
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

#output select ASV count table
asv_tab <- t(gyrb.asv)
dim(asv_tab)
row.names(asv_tab) <- sub(">", "", asv_headers)
select_asv_tab <- asv_tab[row.names(asv_tab) %in% tax_table_ASVs, ]
print('Bacteroidales ASV table')
dim(select_asv_tab)
write.table(select_asv_tab, file.path(outpath,"ASVs_filtered_counts.tsv"),
            sep="\t", quote=F, col.names=NA)

#verify sample names in metadata
all_samples_metadata_file = "metadata/metadata_gyrb_amp_meta.txt"
all_samples_metadata <- read.csv(all_samples_metadata_file,sep='\t')

setdiff(colnames(asv_tab),all_samples_metadata$X.SampleID)
metadata = all_samples_metadata %>%  #filter metadata based on whether samples are in OTU table
  filter(X.SampleID %in% colnames(asv_tab)) %>%
  column_to_rownames('X.SampleID')
write.table(metadata,file = "metadata/metadata_gyrb_amp_meta_passing_samples.txt",sep='\t',quote=FALSE)

#generate phylogeny
paste0('./scripts/gyrb/processing/merge_amp_metagenomic_data/run_phylogeny.sh')

tree <- read_tree(file.path(outpath,"phylogeny/ASVs_filtered_ref_full.fasta.aln.full.tree"))

outgroup_taxa <- c("GCA_001780825.1p__Gemmatimonadotao__GWA2-58-10f__GWA2-58-10s__GWA2-58-10_sp001780825",
                   "GCA_002348465.1p__Gemmatimonadotao__SG8-23f__UBA6960s__UBA2589_sp002348465",
                   "GCA_002686955.1p__Gemmatimonadota_Ao__GCA-2686955f__GCA-2686955s__GCA-2686955_sp002686955")
outgroup_MRCA <- findMRCA(tree,outgroup_taxa)
tree_rooted <- root(tree,node=outgroup_MRCA)
write.tree(tree_rooted,file.path(outpath,"phylogeny/ASVs_filtered_ref_full.fasta.aln.full.rooted.tree"))

tree <- read_tree(file.path(outpath,"phylogeny/ASVs_filtered_ref_full.fasta.aln.amp.tree"))

outgroup_taxa <- c("GCA_001780825.1p__Gemmatimonadotao__GWA2-58-10f__GWA2-58-10s__GWA2-58-10_sp001780825",
                   "GCA_002348465.1p__Gemmatimonadotao__SG8-23f__UBA6960s__UBA2589_sp002348465",
                   "GCA_002686955.1p__Gemmatimonadota_Ao__GCA-2686955f__GCA-2686955s__GCA-2686955_sp002686955")
outgroup_MRCA <- findMRCA(tree,outgroup_taxa)
outgroup_MRCA
tree_rooted <- root(tree,node=outgroup_MRCA)
write.tree(tree_rooted,file.path(outpath,"phylogeny/ASVs_filtered_ref_full.fasta.aln.amp.rooted.tree"))


