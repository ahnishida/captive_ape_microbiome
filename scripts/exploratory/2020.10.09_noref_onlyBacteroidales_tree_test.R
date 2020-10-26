library(dada2); packageVersion("dada2")
library(phyloseq)
library(DECIPHER)
library(tidyverse)
library(phytools)
library(seqinr)
library(Biostrings)
#library("genefilter")
setwd('/Volumes/AHN/captive_ape_microbiome')

#this merges amplicon datasets and filters out reads that aren;t assign. 
#filtered ASVs are then merged with the gtdbtk reference fasta and used to blastn query metagenomic datasets

#load assign taxonomy script
source('scripts/processing/gyrb/idTaxa.R')

#define inputs/outputs
gyrb_ref_faa = 'ref_seqs/gtdbtk_gyrb.faa'
outpath='results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2'
dir.create(outpath, recursive=TRUE)

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
  filter(length>(83*.90)&length<(83*1.10))  %>%
  filter(pident > 80)  %>%
  filter(!ASV %in% ASVs_w_stop_codons)
  
hist(ASVs_filtered_blast$sstart)
hist(ASVs_filtered_blast$send)
ASVs_filtered_blast <-ASVs_filtered_blast %>% 
  filter(sstart<150) %>% 
  filter(send<220)
hist(ASVs_filtered_blast$sstart)
hist(ASVs_filtered_blast$send)

ASVs_filtered_blast_fasta <- asv_fasta[names(asv_fasta) %in% ASVs_filtered_blast$ASV]
Biostrings::writeXStringSet(ASVs_filtered_blast_fasta ,file=file.path(outpath,"ASVs_filtered_blast.fasta")) 

#assign taxonomy
TAX_table <- assign_taxonomy_w_idTAXA(file.path(outpath,"ASVs_all.faa"),gyrb_ref_faa)
TAX_table %>% group_by(Order) %>% tally()   

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

#test refs seqs and alignment softward to see how variation in phylogeny influences results
#get Bacteroidales ref seqs to include with phylogeny
gyrb_ref_fasta = 'ref_seqs/gtdbtk_gyrb.fasta'
gyrb_ref <- readDNAStringSet(gyrb_ref_fasta,format = 'fasta')
str_detect(TREE$tip.label,outgroup)]


alignTrans_ref_ASV_fastas <- function(ASV_fasta,ref_fasta){
  asv_fasta_filt <- Biostrings::readDNAStringSet(file.path(outpath,folder,"ASVs_filtered.fasta"))
  ref <- Biostrings::readDNAStringSet(file.path(outpath,folder,"ref.fasta"))
  asv_ref <- c(ref,DNAStringSet(asv_fasta_filt,start=2))
  Biostrings::writeXStringSet(asv_ref,file=file.path(outpath,folder,"ASVs_filtered_ref.fasta")) 
  asv_ref_aln <- AlignTranslation(asv_ref)  
  Biostrings::writeXStringSet(asv_ref_aln,file=file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta")) 
}

trimAlign <- function(aln_in,aln_out){
   #trim alignment to start and end of ASVs
  aln <- Biostrings::readDNAStringSet(aln_in)
  ASV_1 <- DNAString(as.character(aln['ASV_1']))
  ASV_seq_no_gaps <- as.character(RemoveGaps(aln['ASV_1'],removeGaps = "all"))
  first10nt <- stringr::str_sub(ASV_seq_no_gaps,1,10)
  last10nt <- stringr::str_sub(ASV_seq_no_gaps,-10,-1)
  s <- start(matchPattern(first10nt, ASV_1))
  e <-end(matchPattern(last10nt, ASV_1))
  alnTrim <- DNAStringSet(aln,start=s,end=e)
  seqlengths = width(RemoveGaps(alnTrim,
                                removeGaps = "all",
                                processors = 1))
  alnTrimFilt <- alnTrim[seqlengths > 250*.95]
  alnTrimFilt <- RemoveGaps(alnTrimFilt,
                            removeGaps = "common")
  Biostrings::writeXStringSet(alnTrimFilt,file=aln_out)
}

rootTree <- function(in_tree,rooted_tree,outgroup) {
  TREE <- ape::read.tree(in_tree)
  outgroup_taxa = TREE$tip.label[str_detect(TREE$tip.label,outgroup)]
  mrca <- findMRCA(TREE,outgroup_taxa)
  TREE_rooted <- reroot(TREE,mrca)
  write.tree(TREE_rooted,"tmp.tree")
  system(paste0("sed 's/Root;/;/g' tmp.tree > ",rooted_tree))
  system('rm tmp.tree')
  }


#output ref db fasta tests
gyrb_ref_faa = 'ref_seqs/gtdbtk_gyrb.faa'
gyrb_ref_fasta = 'ref_seqs/gtdbtk_gyrb.fasta'
gyrb_ref <- readDNAStringSet(gyrb_ref_fasta,format = 'fasta')

#generate dataframe with taxnomic info
s <- strsplit(names(gyrb_ref), ";")
gyrb_ref_tax <- data.frame(t(data.frame(s)))
colnames(gyrb_ref_tax) <- c('genome','Phylum','Class','Order','Family','Genus','Species')

#format ref taxa names 
names(gyrb_ref) <- gsub(';','_',names(gyrb_ref))
names(gyrb_ref) <- gsub(' ','__',names(gyrb_ref))
#add fullnames to tax df
gyrb_ref_tax$fullname <- names(gyrb_ref)

gyrb_ref_tax %>% filter(Phylum == "p__Bacteroidota") %>% group_by(Class) %>% tally()
gyrb_ref_tax %>% filter(Class == "c__Bacteroidia") %>% group_by(Order) %>% tally()
gyrb_ref_tax %>% filter(Order == "o__Bacteroidales") %>% group_by(Family) %>% tally()

BacteroidalesOnly_tax <- gyrb_ref_tax %>% filter(Order == "o__Bacteroidales")
BacteroidalesOnly_fasta <- gyrb_ref[names(gyrb_ref) %in% BacteroidalesOnly_tax$fullname]

############
folder <- 'BacteroidalesOnly_AlignTrans'
dir.create(file.path(outpath,'BacteroidalesOnly_AlignTrans'), recursive=TRUE)

#write ref
Biostrings::writeXStringSet(BacteroidalesOnly_fasta,file.path(outpath,folder,'ref.fasta'))

#copy fastas
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
#align ref fasta and ASV fasta
alignTrans(folder)
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('f__Lentimicrobiaceae','f__UBA7960'))
##########
#########
folder <- 'noRef_AlignTrans'
dir.create(file.path(outpath,folder), recursive=TRUE)

#copy fastas
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))

asv <- readDNAStringSet(file.path(outpath,folder,"ASVs_filtered.fasta"),format = 'fasta')
subset <- names(readDNAStringSet(file.path(outpath,"BacteroidalesOnly_AlignTrans","ASVs_filtered_ref_aln_trim.fasta")))
asv_subset <- asv[names(asv) %in% subset]
removed_asvs <- asv[!names(asv) %in% subset]
asv_align <- AlignTranslation(DNAStringSet(asv,start = 2))
Biostrings::writeXStringSet(asv_align,file.path(outpath,folder,'ASVs_filtered_aln.fasta'))
trimAlign(file.path(outpath,folder,'ASVs_filtered_aln.fasta'),
          file.path(outpath,folder,"ASVs_filtered_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_aln_trim.tre")))

TAX_table <- read.table(file.path(outpath,'physeq_Bacteroidales_taxonomy.txt'),header=T)
TAX_table %>% filter(ASV %in% names(asv_align)) %>% group_by(Family) %>% tally()
outgroup_F082 <- TAX_table %>% filter(Family == 'f__F082') 
rootTree(file.path(outpath,folder,"ASVs_filtered_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_aln_trim_rooted.tre"),
         outgroup = outgroup_F082$ASV)
system(paste('cp',file.path(outpath,folder,"ASVs_filtered_aln_trim_rooted.tre"),file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),sep=' '))


