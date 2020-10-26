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

#test refs seqs and alignment softward to see how variation in phylogeny influences results

alignTrans <- function(folder){
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
  first10nt <- stringr::str_sub(ASV_seq_no_gaps,1,8)
  last10nt <- stringr::str_sub(ASV_seq_no_gaps,-8,-1)
  s <- start(matchPattern(first10nt, ASV_1))
  e <-end(matchPattern(last10nt, ASV_1))
  alnTrim <- DNAStringSet(aln,start=s,end=e)
  seqlengths = width(RemoveGaps(alnTrim,
                                removeGaps = "all",
                                processors = 1))
  alnTrimFilt <- alnTrim[seqlengths > 250*.90]
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


folder <- 'GemRef_AlignTrans'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
system(paste0('cp results/gyrb/processing/ref_gyrb_gtdbtk/gyrb_fastas/gtdbtk_gyrb_Bt.fasta ',file.path(outpath,folder,'ref.fasta')))
alignTrans(folder)
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('p__Gemmatimonadota'))
############
folder <- 'BtRef_AlignTrans'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
system(paste0('cp ref_seqs/Bt_out.fasta ',file.path(outpath,folder,'ref.fasta')))
alignTrans(folder)
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('o__NS11-12g'))
############
folder <- 'noRefNSoutgroup_AlignTrans'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
system(paste0('cp ref_seqs/Bt_outgroup_only.fasta ',file.path(outpath,folder,'ref.fasta')))
alignTrans(folder)
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('o__NS11-12g'))
############
############
############
###MAFFT####
folder <- 'GemRef_AlignMafft'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
system(paste0('cp results/gyrb/processing/ref_gyrb_gtdbtk/gyrb_fastas/gtdbtk_gyrb_Bt.fasta ',file.path(outpath,folder,'ref.fasta')))

folder <- 'BtRef_AlignMafft'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
system(paste0('cp ref_seqs/Bt_out.fasta ',file.path(outpath,folder,'ref.fasta')))

folder <- 'noRefNSoutgroup_AlignMafft'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
system(paste0('cp ref_seqs/Bt_outgroup_only.fasta ',file.path(outpath,folder,'ref.fasta')))

#'./scripts/processing/gyrb_old/merge_amp_metagenomic_data/run_phylogeny_10.2.sh '

folder <- 'GemRef_AlignMafft'
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('p__Gemmatimonadota'))

folder <- 'BtRef_AlignMafft'
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('o__NS11-12g'))

folder <- 'noRefNSoutgroup_AlignMafft'
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('o__NS11-12g'))
############


############
folder <- 'BacteroidalesOnly_AlignTrans'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
alignTrans(folder)
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('f__P3','f__UBA7960'))
##########
#########
folder <- 'noRef_AlignTrans'
system(paste0('mkdir ',file.path(outpath,folder)))
system(paste('cp',file.path(outpath,"ASVs_filtered.fasta"),file.path(outpath,folder,"ASVs_filtered.fasta"),sep=' '))
system(paste('cp',file.path(outpath,"ASVs_filtered.faa"),file.path(outpath,folder,"ASVs_filtered.faa"),sep=' '))
system(paste0('touch ',file.path(outpath,folder,'ref.fasta')))
alignTrans(folder)
trimAlign(file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta"),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
rootTree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"),
         file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
         outgroup = c('ASV_8393','ASV_8397'))
