setwd('/Volumes/AHN/captive_ape_microbiome')
library(ape)
library(ggtree)
library(tidytree)
library(tidyverse)
library(phytools)
library(DECIPHER)

#get ASVs of interest and their sequences from old analyses folder
OLD_ASV_table <- read.table("results/gyrb/analyses_old/figures/Figure3_table.txt",sep='\t',header=TRUE)
OLD_CP_ASVs <- OLD_ASV_table %>% filter(CP_pres == 'True')
table(OLD_CP_ASVs$HR_type)

#get ASVs of interest and their sequences from new analyses folder
NEW_ASV_table <- read.table("results/gyrb/analyses/figures/Figure3_table.txt",sep='\t',header=TRUE)
NEW_CP_ASVs <- NEW_ASV_table %>% filter(CP_pres == 'True')
table(NEW_CP_ASVs$HR_type)

old_fasta <- readDNAStringSet(file.path("results/gyrb/inputs_old/ASVs_filtered.fasta"),format = 'fasta') 
new_fasta <- readDNAStringSet("results/gyrb/processing/gyrb_amp_meta_datasets/ASVs_filtered.fasta")
setdiff(new_fasta,old_fasta)
setdiff(old_fasta,new_fasta)

#WHAT HAPPENED TO THE BON GOR CAPTIVE CLADE?
CW_ASVs <- OLD_ASV_table %>% 
  filter(HR_type == 'MX_2_wild_apes'| HR_type == 'HR_wild_chimp') %>% 
  select(ASV, cladeName, HR_sampleTypes)

bon <- OLD_ASV_table %>% filter(cladeName == 'clade_617')
bon_old_fasta <- old_fasta[names(old_fasta) %in% bon$ASV]
bon_new_fasta = intersect(new_fasta,bon_old_fasta)
setdiff(bon_old_fasta,new_fasta) #seqs not in new fasta file
NEW_ASV_table %>% filter(ASV %in% names(bon_new_fasta)) %>% 
  select(ASV, cladeName, HR_sampleTypes,lineage,captive_all,captiveNames)


#WHAT HAPPENED TO THE WILD CHIMP CLADE IN BT2 lineage
ch_ASV_tab <- OLD_CP_ASVs %>% filter(HR_type == 'HR_wild_chimp') %>% 
  select(ASV, cladeName, HR_sampleTypes,captive_all,captiveNames)
asv_fasta <- readDNAStringSet(file.path("results/gyrb/inputs_old/ASVs_filtered.fasta"),format = 'fasta')
ch_ASV_fasta <- asv_fasta[names(asv_fasta) %in%  ch_ASV_tab$ASV]
ch_ASVs_new = intersect(ASVs_raw,ch_ASV_fasta)
NEW_ASV_table %>% filter(ASV %in% names(ch_ASVs_new)) %>% 
  select(ASV, cladeName, HR_sampleTypes,lineage,captive_all,captiveNames)




                         