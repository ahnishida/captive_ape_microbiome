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
outpath='results/gyrb/processing/gyrb_amp_meta_datasets_10.2'
dir.create(outpath, recursive=TRUE)
final_outdir='results/gyrb/inputs_10.2'
dir.create(final_outdir, recursive=TRUE)
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
print("import human_metagenomic_data")
metagenomic_path <- "results/gyrb/processing/metagenomic_samples/final_outputs/gyrb_metagenomic_seqtab.txt"
metagenomic_human_reads <- as.matrix(read.table(metagenomic_path,header=T,
                                                row.names=1, check.names=F, sep="\t")) #Opens ASV table
dim(metagenomic_human_reads)

#merge ASV tables
seqtab <- mergeSequenceTables(gyrb_moeller_wild.seqtab.nochim,
                                gyrb_nishida_captive_wild.seqtab.nochim,
                                gyrb_nishida_projectchimps.seqtab.nochim,
                                metagenomic_human_reads)
saveRDS(seqtab,file.path(outpath,'seqtab.rds'))
seqtab = readRDS(file.path(outpath,'seqtab.rds'))

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(outpath,"ASVs_all.fasta"))
rownames(seqtab)

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)

# metadata
all_samples_metadata_file = "metadata/metadata_gyrb_amp_meta.txt"
metadata <- read.csv(all_samples_metadata_file,sep='\t')
metadata = metadata %>%  #filter metadata based on whether samples are in OTU table
  filter(X.SampleID %in% colnames(asv_tab)) %>%
  column_to_rownames('X.SampleID')
setdiff(colnames(asv_tab),row.names(metadata)) #samples in otu table that aren't in metadata

#filter ASVs
asv_fasta <- readDNAStringSet(file.path(outpath,"ASVs_all.fasta"),format = 'fasta')
asv_faa <- Biostrings::translate(DNAStringSet(asv_fasta,start =2))
Biostrings::writeXStringSet(asv_faa,file=file.path(outpath,"ASVs_all.faa")) 
system(paste('./scripts/processing/gyrb/blastp_filter_ASVs.sh',
      file.path(outpath,"ASVs_all.faa"),
      gyrb_ref_faa, 
      file.path(outpath,"ASVs_all_blastp.txt"),sep=' '))
ASV_blastp_res <- read.table(file.path(outpath,"ASVs_all_blastp.txt"))
colnames(ASV_blastp_res) <- c('ASV','sacc','pident','qlen','length','evalue','bitscore','genome',
                              'salltitles','species','sseq','qseq','sstart','send')

#identity asvs that when translated have stop codons
ASVs_w_stop_codons <- names(asv_faa[str_detect(asv_faa,'\\*')])
length_cutoff = .90
ASVs_filtered_blast <- ASV_blastp_res %>%
  filter(length>(83*length_cutoff))  %>%
  filter(!ASV %in% ASVs_w_stop_codons) 
ASVs_filtered_blast <- ASVs_filtered_blast  %>%
  filter(sstart > .95*median(ASVs_filtered_blast$sstart)) %>%
  filter(sstart < 1.05*median(ASVs_filtered_blast$sstart)) %>%
  filter(send > .95*median(ASVs_filtered_blast$send)) %>%
  filter(send < 1.05*median(ASVs_filtered_blast$send)) 
 
ASVs_filtered_blast_fasta <- asv_fasta[names(asv_fasta) %in% ASVs_filtered_blast$ASV]
Biostrings::writeXStringSet(ASVs_filtered_blast_fasta ,file=file.path(outpath,"ASVs_filtered_blast.fasta")) 

#assign taxonomy
TAX_table <- assign_taxonomy_w_idTAXA(file.path(outpath,"ASVs_all.faa"),gyrb_ref_faa)
table(TAX_table$Order)    

# import into phyloseq
(ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
                sample_data(metadata), 
                tax_table(as.matrix(TAX_table))))
saveRDS(ps,file.path(outpath,'ps.rds'))
ps = readRDS(file.path(outpath,'ps.rds'))

#write physeq
ps <- prune_taxa(ASVs_filtered_blast$ASV,ps)
ps_Bacteroidales <- subset_taxa(ps, Order == 'o__Bacteroidales')
ps_Bacteroidales <- prune_samples(sample_sums(ps_Bacteroidales)>0, ps_Bacteroidales)
saveRDS(ps_Bacteroidales,file.path(outpath,"physeq_Bacteroidales.rds"))

#write fasta
asv_fasta_filt <- asv_fasta[names(asv_fasta) %in% taxa_names(ps_Bacteroidales)]
Biostrings::writeXStringSet(asv_fasta_filt,file=file.path(outpath,"ASVs_filtered.fasta")) 
asv_faa_filt <- asv_faa[names(asv_faa) %in% taxa_names(ps_Bacteroidales)]
Biostrings::writeXStringSet(asv_faa_filt,file=file.path(outpath,"ASVs_filtered.faa")) 

#write to inputs folder
tax_table_tab <- as.data.frame(tax_table(ps_Bacteroidales)) %>% 
  rownames_to_column(var='ASV')
write.table(tax_table_tab,file.path(outpath,'physeq_Bacteroidales_taxonomy.txt'),quote=F,row.names=F,sep='\t')

metadata_tab <- as.data.frame(sample_data(ps_Bacteroidales)) 
metadata_tab$X.SampleID <- rownames(metadata_tab)
write.table(metadata_tab,file.path(outpath,'physeq_metadata_passing_samples.txt'),quote=F,row.names=F,sep='\t')

otu_table <- as.data.frame(otu_table(ps_Bacteroidales)) 
write.table(otu_table,file.path(outpath,'physeq_Bacteroidales_asv_tab.txt'),quote=F,row.names=T,col.names=NA,sep='\t')

#phylogeny
#select ref seqs to include with phylogeny
gyrb_ref_fasta = 'ref_seqs/gtdbtk_gyrb.fasta'
gyrb_ref <- readDNAStringSet(gyrb_ref_fasta,format = 'fasta')
s <- strsplit(names(gyrb_ref), ";")
phylum <- sapply(s, `[`, 2)
class <- sapply(s, `[`, 3)
order <- sapply(s, `[`, 4)
table(class[phylum == "p__Bacteroidota"]) #make c__Chlorobia and c__Ignavibacteria the outgroups
table(order[class == "c__Bacteroidia"]) 
table(order[class == "c__Bacteroidia"]) 
#Bacteroidia orders
Bacteroidia_ref <- gyrb_ref[order %in% c('o__AKYH767-A','o__AKYH767','o__Bacteroidales')]
names(Bacteroidia_ref)
#combine ASVs and ref seqs and align
asv_fasta_filt <- Biostrings::readDNAStringSet(file.path(outpath,"ASVs_filtered.fasta"))
asv_ref <- c(Bacteroidia_ref,DNAStringSet(asv_fasta_filt,start=2))
asv_ref_aln <- AlignTranslation(asv_ref)
Biostrings::writeXStringSet(asv_ref_aln,file=file.path(outpath,"ASVs_filtered_aln.fasta")) 
asv_ref_aln_trim <- DNAStringSet(asv_ref_aln,start=535)
Biostrings::writeXStringSet(asv_ref_aln_trim,file=file.path(outpath,"ASVs_filtered_aln_trim.fasta")) 
#run fasttree
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,"ASVs_filtered_aln_trim.fasta"),
              ' > ',
              file.path(outpath,"ASVs_filtered_aln_trim.tre")))

TREE <- ape::read.tree(file.path(outpath,"ASVs_filtered_aln_trim.tre"))
outgroup_taxa = names(gyrb_ref[order %in% c('o__AKYH767-A','o__AKYH767')])
outgroup_taxa <- strsplit(outgroup_taxa, " ")
outgroup_taxa <- sapply(outgroup_taxa, `[`, 1)
print(outgroup_taxa)
mrca <- findMRCA(TREE,outgroup_taxa)
reroot(TREE,mrca)
write.tree(TREE,(file.path(outpath,"ASVs_filtered_aln_trim_rooted.tre")))
TREE_subset <- keep.tip(TREE,taxa_names(ps_Bacteroidales))
TREE_subset
#add tree to phylogeny
ps_Bacteroidales <- merge_phyloseq(ps_Bacteroidales,TREE_subset)
plot_tree(ps_Bacteroidales, color="Family")
plot_tree(ps_Bacteroidales, color="Order")

write.tree(phy_tree(ps_Bacteroidales),file.path(outpath,'physeq_Bacteroidales.tree'))  
write.tree(phy_tree(ps_Bacteroidales),file.path(final_outdir,'physeq_Bacteroidales.tree'))  

Biostrings::writeXStringSet(asv_fasta_filt,file=file.path(final_outdir,"physeq_Bacteroidales_asv.fasta")) 


