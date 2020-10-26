library(dada2); packageVersion("dada2")
library(phyloseq)
library(DECIPHER)
library(tidyverse)
library(phytools)
library(seqinr)
library(Biostrings)
#library("genefilter")
setwd('/Volumes/AHN/captive_ape_microbiome')

#COMPARE Gemma Ref to other no

setwd('results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2')

BtRef <- Biostrings::readDNAStringSet(file.path("BtRef_AlignTrans","ref.fasta"))
GemRef <- Biostrings::readDNAStringSet(file.path("GemRef_AlignTrans","ref.fasta"))

common_seqs <- BtRef[BtRef %in% GemRef]

#unique to Bt Res
BtRef_unique <- BtRef[!BtRef %in% GemRef]
names(BtRef_unique)

#unique to Bt Res
GemRef_unique <- GemRef[!GemRef %in% BtRef]
names(GemRef_unique)

#output ref db fasta tests
setwd('/Volumes/AHN/captive_ape_microbiome')
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
dir.create(file.path('results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2','BacteroidalesOnly_AlignTrans'), recursive=TRUE)
Biostrings::writeXStringSet(BacteroidalesOnly_fasta,file.path('results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2','BacteroidalesOnly_AlignTrans','ref.fasta'))

Bt_Family5_ref_tax <- gyrb_ref_tax %>% filter(Order == "o__Bacteroidales") %>% 
  group_by(Family) %>% 
  slice(1:5)
table(Bt_Family5_ref_tax$Family) 
Bt_Family5_ref_fasta <- gyrb_ref[names(gyrb_ref) %in% Bt_Family5_ref_tax$fullname]
dir.create(file.path('results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2','Bt_Family5_ref'), recursive=TRUE)
Biostrings::writeXStringSet(Bt_Family5_ref_fasta,file.path('results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2','Bt_Family5_ref','ref.fasta'))

# <- 
TAX_table_file <- file.path('results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2','physeq_Bacteroidales_taxonomy.txt')
TAX_table <- read.table(TAX_table_file,header=TRUE)
table(TAX_table$Family)
TAX_table %>% filter(Family == 'f__F082')
Bt_select_fam <- gyrb_ref[family %in% unique(TAX_table$Family)]
Bt_select_fam 


