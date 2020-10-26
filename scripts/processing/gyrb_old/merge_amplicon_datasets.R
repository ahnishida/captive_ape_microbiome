library(dada2); packageVersion("dada2")
library(phyloseq)
library(DECIPHER)
library(tidyverse)
library(phytools)
library(seqinr)
#library("genefilter")
setwd('/Volumes/AHN/captive_ape_microbiome')

#this merges amplicon datasets and filters out reads that aren;t assign. 
#filtered ASVs are then merged with the gtdbtk reference fasta and used to blastn query metagenomic datasets

#load assign taxonomy script
source('scripts/processing/gyrb/idTaxa.R')

#define inputs/outputs
gyrb_ref_faa = 'ref_seqs/gtdbtk_gyrb.faa'
outpath='results/gyrb/processing/gyrb_amp_datasets'
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
seqtab <- mergeSequenceTables(gyrb_moeller_wild.seqtab.nochim,
                                gyrb_nishida_captive_wild.seqtab.nochim,
                                gyrb_nishida_projectchimps.seqtab.nochim)


# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# correcting sample names so they don't have funny extensions
samples.out <- rownames(seqtab)
samples.out <- sapply(strsplit(basename(samples.out), "_"), `[`, 1)
rownames(seqtab) <- samples.out

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)


# making and writing out a fasta of ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(outpath,"ASVs_all.fasta"))

#filter ASVs
asv_fasta <- readDNAStringSet(file.path(outpath,"ASVs_all.fasta"),format = 'fasta')
asv_faa <- Biostrings::translate(DNAStringSet(asv_fasta,start =2))
Biostrings::writeXStringSet(asv_faa,file=file.path(outpath,"ASVs_all.faa")) 
system(paste('./scripts/processing/gyrb/blastp_filter_ASVs.sh',
      file.path(outpath,"ASVs_all.faa"),
      gyrb_ref_faa, 
      file.path(outpath,"ASVs_all_blastp.txt"),sep=' '))
ASV_blastp_res <- read.table(file.path(outpath,"ASVs_all_blastp.txt"))
colnames(ASV_blastp_res) <- c('qseqid','sseqid','pident','length','qlen','evalue')
hist(ASV_blastp_res$pident)
ASV_filtered <- ASV_blastp_res %>%  filter(pident > 80 & length > qlen*.90 & length < qlen*1.10)
asv_fasta_filt <- asv_fasta[names(asv_fasta) %in% ASV_filtered$qseqid]
Biostrings::writeXStringSet(asv_fasta_filt,file=file.path(outpath,"ASVs_filtered.fasta")) 
asv_faa_filt <- asv_faa[names(asv_faa) %in% ASV_filtered$qseqid]
Biostrings::writeXStringSet(asv_faa_filt,file=file.path(outpath,"ASVs_filtered.faa")) 

#assign taxonomy
TAX_table <- assign_taxonomy_w_idTAXA(file.path(outpath,"ASVs_filtered.faa"),gyrb_ref_faa)
table(TAX_table$Class)    

# metadata
all_samples_metadata_file = "metadata/metadata_gyrb_amp_meta.txt"
metadata <- read.csv(all_samples_metadata_file,sep='\t')
metadata = metadata %>%  #filter metadata based on whether samples are in OTU table
  filter(X.SampleID %in% colnames(asv_tab)) %>%
  column_to_rownames('X.SampleID')
setdiff(colnames(asv_tab),row.names(metadata)) #samples in otu table that aren't in metadata

# import into phyloseq
(ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
                sample_data(metadata), 
                tax_table(as.matrix(TAX_table))))

#phylogeny
#select ref seqs to include with phylogeny
gyrb_ref_fasta = 'ref_seqs/gtdbtk_gyrb.fasta'
gyrb_ref <- readDNAStringSet(gyrb_ref_fasta,format = 'fasta')
s <- strsplit(names(gyrb_ref), ";")
phylum <- sapply(s, `[`, 2)
class <- sapply(s, `[`, 3)
order <- sapply(s, `[`, 4)
table(class[phylum == "p__Bacteroidota"]) #make c__Chlorobia and c__Ignavibacteria the outgroups
table(order[class == "c__Bacteroidia"]) #Bacteroidia orders
Bacteroidia_ref <- gyrb_ref[class %in% c('c__Chlorobia','c__Ignavibacteria','c__Bacteroidia')]

#combine ASVs and ref seqs and align
asv_fasta_filt <- Biostrings::readDNAStringSet(file.path(outpath,"ASVs_filtered.fasta"))
asv_ref <- c(Bacteroidia_ref,DNAStringSet(asv_fasta_filt,start=2))
asv_ref_aln <- AlignTranslation(asv_ref)
Biostrings::writeXStringSet(asv_ref_aln,file=file.path(outpath,"ASVs_filtered_aln.fasta")) 
asv_ref_aln_trim <- DNAStringSet(asv_ref_aln,start=526)
Biostrings::writeXStringSet(asv_ref_aln_trim,file=file.path(outpath,"ASVs_filtered_aln_trim.fasta")) 
#run fasttree
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,"ASVs_filtered_aln_trim.fasta"),
              ' > ',
              file.path(outpath,"ASVs_filtered_aln_trim.tre")))

TREE <- ape::read.tree(file.path(outpath,"ASVs_filtered_aln_trim.tre"))
outgroup_taxa = names(gyrb_ref[class %in% c('c__Chlorobia','c__Ignavibacteria')])
outgroup_taxa <- strsplit(outgroup_taxa, " ")
outgroup_taxa <- sapply(outgroup_taxa, `[`, 1)
print(outgroup_taxa)
mrca <- findMRCA(TREE,outgroup_taxa)
reroot(TREE,mrca)
write.tree(TREE,(file.path(outpath,"ASVs_filtered_aln_trim_rooted.tre")))
TREE_subset <- keep.tip(TREE,taxa_names(ps))

#add tree to phylogeny
ps <- merge_phyloseq(ps,TREE_subset)
plot_tree(ps, color="Family")
plot_tree(ps, color="Order")

#write to inputs folder
saveRDS(ps,file.path(outpath,"physeq.rds"))
tax_table_tab <- as.data.frame(tax_table(ps)) %>% 
  rownames_to_column(var='ASV')
write.table(tax_table_tab,file.path(outpath,'physeq_taxonomy.txt'),quote=F,row.names=F,sep='\t')

metadata_tab <- as.data.frame(sample_data(ps)) 
metadata_tab$X.SampleID <- rownames(metadata_tab)
write.table(metadata_tab,file.path(outpath,'physeq_metadata_passing_samples.txt'),quote=F,row.names=F,sep='\t')

otu_table <- as.data.frame(otu_table(ps)) 
write.table(otu_table,file.path(outpath,'physeq_asv_tab.txt'),quote=F,row.names=T,col.names=NA,sep='\t')

write.tree(phy_tree(ps),file.path(outpath,'physeq.tree'))  

