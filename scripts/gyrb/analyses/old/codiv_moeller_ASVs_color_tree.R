setwd('/Volumes/AHN/captive_ape_microbiome')
library(ape)
library(ggtree)
library(tidytree)
library(tidyverse)
library(ggrepel)
library(phytools)

full_tree <- read.tree('results/gyrb/inputs/ASVs_filtered_ref_full.tree')
codiv_table <- read.table('results/gyrb/analyses/codiv_moeller_ASVs/codiv_clades_ASVs.txt',sep='\t',header=TRUE)
HRclades_fulltree_ASVs <- read.table("results/gyrb/analyses/codiv_moeller_ASVs/HRclades_fulltree_ASVs.txt",sep='\t',header=TRUE)
unique(codiv_table$codiv_clade)
unique(codiv_table$lineage)

#Subset tree and codivASVs to Bt1 lineage
lin = 'Bt1'
lin_table <- codiv_table %>% filter(lineage == 'Bt1')
lin_ASVs = as.vector(lin_table$ASV)
lin_MRCA <- findMRCA(full_tree,lin_ASVs)
lin_tree <- extract.clade(full_tree,lin_MRCA)
ggtree(lin_tree) + geom_nodepoint(aes(subset = label > .50))


#moeller clade labels
lin_clades <- as.vector(unique(lin_table$codiv_clade))
nnode <- c()
for (clade in lin_clades) {
  clade_ASVs <- lin_table %>% filter(codiv_clade == clade)
  clade_MRCA <- findMRCA(lin_tree,as.vector(clade_ASVs$ASV))
  nnode <- c(nnode,clade_MRCA)
}
codiv_clades <- cbind(lin_clades,nnode)
codiv_clades
ggtree(lin_tree) + #xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=43,label='Bt1_gorilla',color='green4') + 
  geom_cladelabel(node=53,label='Bt1_chimp',color='orange') + 
  geom_cladelabel(node=60,label='Bt1_bonobo',color='red') 

#import in HRclades filter to those present in lineage
HRclades_lin_table <- HRclades_fulltree_ASVs %>% filter(ASV %in% lin_ASVs)
HRclades_lin <- HRclades_lin_table %>% select(cladeName,HR_clade,HR_sampleNum) %>% unique()
(HRcladeNames_lin = as.vector(unique(HRclades_lin$cladeName)))
cls <- list()
for (clade in HRcladeNames_lin) {
  clade_table <- HRclades_lin_table %>% filter(cladeName == clade)
  clade_ASVs <- as.character(clade_table$ASV)
  cls[[length(cls)+1]] <- clade_ASVs
}
#add HR clades to tree
lin_tree_cls <- groupOTU(lin_tree, cls)
color_vec <- as.vector(recode(HRclades_lin$HR_clade, wild_gorilla = "green4", 
                              wild_chimp = "orange",
                              wild_bonobo = "red", 
                              wild_human = "blue"))
(color_vec <- c("black",color_vec))

ggtree(lin_tree_cls, aes(color=group)) +
  scale_color_manual(values=color_vec)

ggtree(lin_tree_cls, aes(color=group)) + xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=43,label='Bt1_gorilla',color='green4') + 
  geom_cladelabel(node=53,label='Bt1_chimp',color='orange') + 
  geom_cladelabel(node=60,label='Bt1_bonobo',color='red') +
  scale_color_manual(values=color_vec)
ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt1_tree.pdf')


#Subset tree and codivASVs to Bt2 lineage
lin_table <- codiv_table %>% filter(lineage == 'Bt2')
lin_ASVs = as.vector(lin_table$ASV)
lin_MRCA <- findMRCA(full_tree,lin_ASVs)
lin_tree <- extract.clade(full_tree,lin_MRCA)
ggtree(lin_tree) + geom_nodepoint(aes(subset = label > .50))


#moeller clade labels
lin_codiv_table<- lin_table %>% filter(codiv_clade != '')
(lin_clades <- as.vector(unique(lin_codiv_table$codiv_clade)))

nnode <- c()
for (clade in lin_clades) {
  clade_ASVs <- lin_table %>% filter(codiv_clade == clade)
  clade_MRCA <- findMRCA(lin_tree,as.vector(clade_ASVs$ASV))
  nnode <- c(nnode,clade_MRCA)
}
codiv_clades <- cbind(lin_clades,nnode)
codiv_clades
ggtree(lin_tree) + xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .50),size = .75) +
  geom_cladelabel(node=402,label='Bt2_clade2_bonobo',color='red') + 
  geom_cladelabel(node=435,label='Bt2_clade2_chimp',color='orange') + 
  geom_cladelabel(node=511,label='Bt2_clade2_bonobo',color='red') +
  geom_cladelabel(node=600,label='Bt2_clade1_chimp',color='orange') + 
  geom_cladelabel(node=626,label='Bt2_clade2_gorilla',color='green4') 

#import in HRclades filter to those present in lineage
HRclades_lin_table <- HRclades_fulltree_ASVs %>% filter(ASV %in% lin_ASVs)
HRclades_lin <- HRclades_lin_table %>% select(cladeName,HR_clade,HR_sampleNum) %>% unique()
(HRcladeNames_lin = as.vector(unique(HRclades_lin$cladeName)))
cls <- list()
for (clade in HRcladeNames_lin) {
  clade_table <- HRclades_lin_table %>% filter(cladeName == clade)
  clade_ASVs <- as.character(clade_table$ASV)
  cls[[length(cls)+1]] <- clade_ASVs
}
#add HR clades to tree
lin_tree_cls <- groupOTU(lin_tree, cls)
color_vec <- as.vector(recode(HRclades_lin$HR_clade, wild_gorilla = "green4", 
                              wild_chimp = "orange",
                              wild_bonobo = "red", 
                              human = "blue"))
(color_vec <- c("black",color_vec))

ggtree(lin_tree_cls, aes(color=group)) +
  scale_color_manual(values=color_vec)

ggtree(lin_tree_cls, aes(color=group)) + xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .70),size=.75) +
  geom_cladelabel(node=402,label='Bt2_clade2_bonobo',color='red') + 
  geom_cladelabel(node=435,label='Bt2_clade2_chimp',color='orange') + 
  geom_cladelabel(node=511,label='Bt2_clade2_bonobo',color='red') +
  geom_cladelabel(node=600,label='Bt2_clade1_chimp',color='orange') + 
  geom_cladelabel(node=626,label='Bt2_clade2_gorilla',color='green4') +
  scale_color_manual(values=color_vec)
ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt2_tree.pdf')

#Subset tree and codivASVs to Bt3 lineage
lin_table <- codiv_table %>% filter(lineage == 'Bt3')
lin_ASVs = as.vector(lin_table$ASV)
lin_MRCA <- findMRCA(full_tree,lin_ASVs)
lin_tree <- extract.clade(full_tree,lin_MRCA)
ggtree(lin_tree) + geom_nodepoint(aes(subset = label > .50))


#moeller clade labels
lin_codiv_table<- lin_table %>% filter(codiv_clade == 'Bt3_clade1_human'|codiv_clade == 'Bt3_clade1_chimp'|codiv_clade == 'Bt3_clade1_bonobo')
(lin_clades <- as.vector(unique(lin_codiv_table$codiv_clade)))

nnode <- c()
for (clade in lin_clades) {
  clade_ASVs <- lin_table %>% filter(codiv_clade == clade)
  clade_MRCA <- findMRCA(lin_tree,as.vector(clade_ASVs$ASV))
  nnode <- c(nnode,clade_MRCA)
}
codiv_clades <- cbind(lin_clades,nnode)
codiv_clades
ggtree(lin_tree) + 
  geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=7999,label='Bt3_clade1_human',color='blue') +
  geom_cladelabel(node=12203,label='Bt3_clade1_chimp',color='orange') 
 # geom_cladelabel(node=7996,label='Bt3_clade1_bonobo',color='red') 

#import in HRclades filter to those present in lineage
HRclades_lin_table <- HRclades_fulltree_ASVs %>% filter(ASV %in% lin_ASVs)
HRclades_lin <- HRclades_lin_table %>% select(cladeName,HR_clade,HR_sampleNum) %>% unique()
(HRcladeNames_lin = as.vector(unique(HRclades_lin$cladeName)))
cls <- list()
for (clade in HRcladeNames_lin) {
  clade_table <- HRclades_lin_table %>% filter(cladeName == clade)
  clade_ASVs <- as.character(clade_table$ASV)
  cls[[length(cls)+1]] <- clade_ASVs
}
#add HR clades to tree
lin_tree_cls <- groupOTU(lin_tree, cls)
color_vec <- as.vector(recode(HRclades_lin$HR_clade, wild_gorilla = "green4", 
                              wild_chimp = "orange",
                              wild_bonobo = "red", 
                              human = "blue"))
(color_vec <- c("black",color_vec))

ggtree(lin_tree_cls, aes(color=group)) +
  scale_color_manual(values=color_vec)

ggtree(lin_tree_cls, aes(color=group)) + 
  geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=7999,label='Bt3_clade1_human',color='blue') +
  geom_cladelabel(node=12203,label='Bt3_clade1_chimp',color='orange') +
  geom_cladelabel(node=7996,label='Bt3_clade1_bonobo',color='red',offset = .5) +
  scale_color_manual(values=color_vec)
ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt3_tree.pdf',dix)
