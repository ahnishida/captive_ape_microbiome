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

#add captive ASVs to tree
captive_ASVs <- read.table('results/gyrb/analyses/codiv_moeller_ASVs/Table_sampletypes_site_captive_ASV_count.txt',sep='\t',header=T)
captive_ASVs <- captive_ASVs %>% column_to_rownames(var='ASV') 
captive_ASVs$ASVsums <- rowSums(captive_ASVs)
captive_ASVs <- captive_ASVs %>% 
  filter(ASVsums > 0) %>% 
  rownames_to_column(var='ASV') 


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
color_vec <- as.vector(recode(HRclades_lin$HR_clade, wild_gorilla = "darkgreen", 
                              wild_chimp = "darkorange2",
                              wild_bonobo = "red2", 
                              human = "dodgerblue3"))
(color_vec <- c("black",color_vec))

ggtree(lin_tree_cls, aes(color=group)) +
  scale_color_manual(values=color_vec)

ggtree(lin_tree_cls, aes(color=group)) + xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=43,label='Bt1 gorilla',color='darkgreen') + 
  geom_cladelabel(node=53,label='Bt1 chimp',color='darkorange2') + 
  geom_cladelabel(node=60,label='Bt1 bonobo',color='red2') +
  scale_color_manual(values=color_vec)

#any captive ape ASVs in linage
captive_ASVs %>%
  filter(ASV %in% lin_ASVs) %>%
  nrow()

ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt1_tree.pdf')

#testing how to color a tree
tree <- ggtree(lin_tree, aes(color=HR_clade)) + xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .50)) 
dd <- HRclades_lin_table %>% 
  select(ASV,HR_clade) %>%
  filter(ASV %in% c('ASV_4444','ASV_4172'))
tree %<+% dd + geom_tippoint(aes(color=HR_clade), size=3) 
tree %<+% dd + geom_tippoint(aes(subset=label %in% c('ASV_4444','ASV_4172'),label=label, size=3))


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
  geom_cladelabel(node=402,label='Bt2_bonobo_clade2',color='red') + 
  geom_cladelabel(node=435,label='Bt2_chimp_clade2',color='orange') + 
  geom_cladelabel(node=511,label='Bt2_bonobo_clade1',color='red') +
  geom_cladelabel(node=600,label='Bt2_chimp_clade1',color='orange') + 
  geom_cladelabel(node=626,label='Bt2_gorilla_clade1',color='green4') 

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
color_vec <- as.vector(recode(HRclades_lin$HR_clade, wild_gorilla = "darkgreen", 
                              wild_chimp = "darkorange2",
                              wild_bonobo = "red2", 
                              human = "dodgerblue3"))
(color_vec <- c("black",color_vec))

ggtree(lin_tree_cls, aes(color=group)) +
  scale_color_manual(values=color_vec)

ggtree(lin_tree_cls, aes(color=group)) + xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .70),size=.75) +
  geom_cladelabel(node=402,label='Bt2_clade2_bonobo',color='red2') + 
  geom_cladelabel(node=435,label='Bt2_clade2_chimp',color='darkorange2') + 
  geom_cladelabel(node=511,label='Bt2_clade2_bonobo',color='red2') +
  geom_cladelabel(node=600,label='Bt2_clade1_chimp',color='darkorange2') + 
  geom_cladelabel(node=626,label='Bt2_clade2_gorilla',color='green4') +
  scale_color_manual(values=color_vec)

#add captive ASVs to tree
captive_ASVs %>%
  filter(ASV %in% lin_ASVs) %>%
  nrow()
captive_ASVs_lin <- captive_ASVs %>%
  filter(ASV %in% lin_ASVs)
captive_ASVs_lin <- HRclades_lin_table %>% 
  select(ASV,cladeName,HR_clade) %>% 
  right_join(captive_ASVs_lin,on='ASV')
captive_ASVs_lin
ggtree(lin_tree_cls) %<+% captive_ASVs_lin + 
  geom_tippoint(aes(subset=label %in% captive_ASVs_lin$ASV,label=label,color='red'))


ggtree(lin_tree_cls, aes(color=group)) + xlim(NA, .3) +
  #geom_nodepoint(aes(subset = label > .70),size=.75) +
  geom_cladelabel(node=402,label='Bt2_clade2_bonobo',color='red') + 
  geom_cladelabel(node=435,label='Bt2_clade2_chimp',color='orange') + 
  geom_cladelabel(node=511,label='Bt2_clade2_bonobo',color='red') +
  geom_cladelabel(node=600,label='Bt2_clade1_chimp',color='orange') + 
  geom_cladelabel(node=626,label='Bt2_clade2_gorilla',color='green4') +
  geom_cladelabel(node=373,label='captive_ape_ASVs',color='magenta') +
  scale_color_manual(values=color_vec)+
  geom_tippoint(aes(subset=label %in% captive_ASVs_lin$ASV,label=label))
  
captive_ASVs_lin <- captive_ASVs %>%
  filter(ASV %in% lin_ASVs)
captive_ASVs_lin <- HRclades_lin_table %>% 
  select(ASV,cladeName,HR_clade) %>% 
  right_join(captive_ASVs_lin,on='ASV') 
captive_ASVs_lin
ggtree(lin_tree_cls, aes(color=group)) + xlim(NA, .4) +
  #geom_nodepoint(aes(subset = label > .70),size=.75) +
  geom_tiplab(aes(subset=label %in% captive_ASVs_lin$ASV,
                  label=label), align = TRUE, color='black') +
  geom_cladelabel(node=402,label='Bt2 bonobo clade2',color='red') + 
  geom_cladelabel(node=435,label='Bt2 chimp clade2',color='orange') + 
  geom_cladelabel(node=511,label='Bt2 bonobo clade1',color='red') +
  geom_cladelabel(node=600,label='Bt2 chimp clade1',color='orange') + 
  geom_cladelabel(node=626,label='Bt2_gorilla clade1',color='green4') +
#  geom_cladelabel(node=373,label='captive_ape_ASVs',color='magenta') +
  scale_color_manual(values=color_vec) +
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_gorilla_COLZ>0],
    label=label), size=1.5,x=.25,color='darkolivegreen3',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_chimp_PC>0],
    label=label), size=1.5,x=.28,color='tan1',shape=17) 
  

ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt2_tree.pdf')

#Subset tree and codivASVs to Bt3 lineage
lin_table <- codiv_table %>% filter(lineage == 'Bt3')
lin_ASVs = as.vector(lin_table$ASV)
lin_MRCA <- findMRCA(full_tree,lin_ASVs)
lin_tree <- extract.clade(full_tree,lin_MRCA)
ggtree(lin_tree) + geom_nodepoint(aes(subset = label > .50))

#moeller clade labels
lin_codiv_table <- lin_table %>% 
  filter(codiv_clade != "")
(lin_clades <- as.vector(unique(lin_codiv_table$codiv_clade)))

nnode <- c()
for (clade in lin_clades) {
  clade_ASVs <- lin_codiv_table %>% filter(codiv_clade == clade)
  clade_MRCA <- findMRCA(lin_tree,as.vector(clade_ASVs$ASV))
  nnode <- c(nnode,clade_MRCA)
}
codiv_clades <- cbind(lin_clades,nnode)
ggtree(lin_tree) + 
  #geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=625,label='Bt3_clade1_human',color='blue') +
  geom_cladelabel(node=1057,label='Bt3_clade1_chimp',color='orange') +
  geom_cladelabel(node=1077,label='Bt3_clade1_bonobo',color='red') 

#import in HRclades filter to those present in lineage
HRclades_lin_table <- HRclades_fulltree_ASVs %>% filter(ASV %in% lin_ASVs)
HRclades_lin <- HRclades_lin_table %>% select(cladeName,HR_clade,HR_sampleNum) %>%
#  filter((HR_clade=='human'&HR_sampleNum>100)|HR_clade!='human') %>% 
  unique()

(HRcladeNames_lin = as.vector(unique(HRclades_lin$cladeName)))
cls <- list()
for (clade in HRcladeNames_lin) {
  clade_table <- HRclades_lin_table %>% filter(cladeName == clade)
  clade_ASVs <- as.character(clade_table$ASV)
  cls[[length(cls)+1]] <- clade_ASVs
}
#add HR clades to tree
lin_tree_cls <- groupOTU(lin_tree, cls)
color_vec <- as.vector(recode(HRclades_lin$HR_clade, wild_gorilla = "darkgreen", 
                              wild_chimp = "darkorange2",
                              wild_bonobo = "red2", 
                              human = "dodgerblue3"))
(color_vec <- c("black",color_vec))

ggtree(lin_tree_cls, aes(color=group)) +
  scale_color_manual(values=color_vec)

ggtree(lin_tree_cls, aes(color=group)) + 
  geom_nodepoint(aes(subset = label > .50),size=.75) +
  geom_cladelabel(node=625,label='Bt3_clade1_human',color='dodgerblue3') +
  geom_cladelabel(node=1057,label='Bt3_clade1_chimp',color='darkorange2') +
  geom_cladelabel(node=1077,label='Bt3_clade1_bonobo',color='red2') +
  scale_color_manual(values=color_vec)

#add captive ASVs to tree
captive_ASVs_lin <- captive_ASVs %>%
  filter(ASV %in% lin_ASVs)
captive_ASVs_lin <- HRclades_lin_table %>% 
  select(ASV,cladeName,HR_clade) %>% 
  right_join(captive_ASVs_lin,on='ASV') 
tree_color_scale
ggtree(lin_tree_cls, aes(color=group)) + xlim(NA, 1) +
  geom_nodepoint(aes(subset = label > .5),size=.75) +
  scale_color_manual(values=color_vec) +
  geom_tiplab(aes(subset=label %in% captive_ASVs_lin$ASV,
                  label=label), align = TRUE, color='black') +
  geom_cladelabel(node=625,label='Bt3 human',color='dodgerblue3',offset=.27) +
  geom_cladelabel(node=1057,label='Bt3 chimp',color='darkorange2',offset=.23) +
  geom_cladelabel(node=1077,label='Bt3 bonobo',color='red2',offset=.23) +
  geom_tippoint(aes(
    subset=label == 'ASV_934',
    label=label), size=3,color='tan1',x=.42,shape=17)  + 
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_bonobo_COLZ>0],
    label=label), size=3,x=.32,color='indianred2') + 
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_gorilla_COLZ>0],
    label=label), size=3,x=.345,color='darkolivegreen3') + 
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_orangutan_COLZ>0],
    label=label), size=3,x=.37,color='plum3') +
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_chimp_HOUZ>0],
    label=label), size=3,x=.395,color='tan1',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_gorilla_HOUZ>0],
    label=label), size=3,x=.42,color='darkolivegreen3',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% captive_ASVs_lin$ASV[captive_ASVs_lin$captive_orangutan_HOUZ>0],
    label=label), size=3,x=.455,color='plum3',shape=15) 
moon_color_scale

ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt3_tree.pdf')
colnames(captive_HR)
captive_HR <- HRclades_fulltree_ASVs %>% select(-sampleNames,-ASVsTax) %>% filter(ASV %in% captive_ASVs$ASV) 
table(captive_HR$HR_clade)
captive_HR %>% filter(HR_clade == 'wild_chimp') 

length(unique(HRclades_fulltree_ASVs$cladeName))
