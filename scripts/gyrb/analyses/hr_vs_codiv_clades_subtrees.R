setwd('/Volumes/AHN/captive_ape_microbiome')
library(ape)
library(ggtree)
library(tidytree)
library(tidyverse)
library(ggrepel)
library(phytools)


full_tree <- read.tree('results/gyrb/inputs/ASVs_filtered_ref_full.tree')
ASV_table <- read.table("results/gyrb/analyses/codiv_moeller_ASVs/full_tree_ASV_table.txt",sep='\t',header=TRUE)
Bt2_inset_tree <- read.tree('results/gyrb/inputs/moeller_codiv_lin_Bt2.tree')
Bt3_inset_tree <- read.tree('results/gyrb/inputs/moeller_codiv_lin_Bt3.tree')
inset_table <- read.table("results/gyrb/inputs/moeller_codiv_HRclades.txt",sep='\t',header=TRUE)

#Subset tree and codivASVs to Bt1 lineage
Bt1_codiv_table <- ASV_table %>% filter(lineage == 'Bt1')
Bt1_codiv_ASVs = as.vector(Bt1_codiv_table$ASV)
Bt1_MRCA <- findMRCA(full_tree,Bt1_codiv_ASVs)
Bt1_tree <- extract.clade(full_tree,Bt1_MRCA)
Bt1_all_ASVs <- Bt1_tree$tip.label
Bt1_ASV_table <- ASV_table %>% filter(ASV %in% Bt1_all_ASVs)
ggtree(Bt1_tree) + geom_nodepoint(aes(subset = label > .50))

#moeller clade labels
Bt1_codiv_clades <- as.vector(unique(Bt1_ASV_table$codiv_clade))

get_codiv_table <- function(tree,clades){
  nnode <- c()
  for (clade in clades){
    clade_ASVs <- ASV_table %>% filter(codiv_clade == clade)
    clade_MRCA <- findMRCA(tree,as.vector(clade_ASVs$ASV))
    nnode <- c(nnode,clade_MRCA)
  }
  codiv_clades <- cbind(nnode,clades)
  return(codiv_clades)
}
(Bt1_codiv_clades <- get_codiv_table(Bt1_tree,Bt1_codiv_clades))

ggtree(Bt1_tree) + #xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=43,label='Bt1_gorilla',color='green4') + 
  geom_cladelabel(node=53,label='Bt1_chimp',color='orange') + 
  geom_cladelabel(node=60,label='Bt1_bonobo',color='red') 

#import in HRclades filter to those present in lineage
add_HRclades_to_tree <- function(lineage_tree,table){
  (taxa<-lineage_tree$tip.label) 
  (lineage_table <- table %>% 
      filter(ASV %in% taxa & cladeName != "") %>% 
      group_by(cladeName) %>% 
      sample_n(1))
  (cladeNames = as.vector(lineage_table$cladeName))
  (cladeHostSp = as.vector(lineage_table$HR_clade))
  cls <- list()
  for (clade in cladeNames) {
    print(clade)
    (clade_table <- table %>% filter(cladeName == clade))
    (clade_ASVs <- as.character(clade_table$ASV))
    cls[[length(cls)+1]] <- clade_ASVs
  }
  cls
  tree_cls <- groupOTU(lineage_tree, cls)
  
  color_vec <- as.vector(recode(cladeHostSp, wild_gorilla = "darkgreen", 
                                wild_chimp = "darkorange2",
                                wild_bonobo = "red2", 
                                human = "dodgerblue3"))
  (color_vec <- c("black",color_vec))
  return(c(tree_cls,color_vec))
}
res <- add_HRclades_to_tree(Bt1_tree,ASV_table)
tree_cls <- res[[1]]
color_vec <- res[[2]]

ggtree(tree_cls, aes(color=group)) +
  scale_color_manual(values=color_vec)

ggtree(tree_cls, aes(color=group)) + xlim(NA, .3) +
  geom_nodepoint(aes(subset = label > .50)) +
  geom_cladelabel(node=43,label='Bt1 gorilla',color='darkgreen') + 
  geom_cladelabel(node=53,label='Bt1 chimp',color='darkorange2') + 
  geom_cladelabel(node=60,label='Bt1 bonobo',color='red2') +
  scale_color_manual(values=color_vec)

#any captive ape ASVs in linage


ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt1_tree.pdf')

#Subset tree and codivASVs to Bt2 lineage
Bt2_codiv_table <- ASV_table %>% filter(lineage == 'Bt2')
Bt2_codiv_ASVs = as.vector(Bt2_codiv_table$ASV)
Bt2_MRCA <- findMRCA(full_tree,Bt2_codiv_ASVs)
Bt2_tree <- extract.clade(full_tree,Bt2_MRCA)
Bt2_all_ASVs <- Bt2_tree$tip.label
Bt2_ASV_table <- ASV_table %>% filter(ASV %in% Bt2_all_ASVs)

Bt2_codiv_clades <- as.vector(unique(Bt2_codiv_table$codiv_clade))
(Bt2_codiv_clades <- get_codiv_table(Bt2_tree,Bt2_codiv_clades))
res <- add_HRclades_to_tree(Bt2_tree,ASV_table)
Bt2_tree_cls <- res[[1]]
Bt2_color_vec <- res[[2]]
(Bt2_captive_ASVs <- ASV_table %>% 
  filter(ASV %in% Bt2_tree$tip.label) %>% 
  filter(captive_all>0))
Bt2_captive_ASVs %>%
  select(ASV,captive_chimp_HOUZ:captive_chimp_PC) %>%
  gather(key = "host_cat", value="sample_count",-ASV) %>%
  filter(sample_count>0)

ggtree(Bt2_tree_cls, aes(color=group)) +
  scale_color_manual(values=Bt2_color_vec)

ggtree(Bt2_tree_cls, aes(color=group)) + xlim(NA, .4) +
  scale_color_manual(values=Bt2_color_vec) + #HR clades
  geom_treescale(x=0,y=0) +
  #geom_nodepoint(aes(subset = label > .70),size=.75) + #bootstrap
  #moeller codiv clades
  geom_cladelabel(node=402,label='Bt2_clade2_bonobo',color='red2') + 
  geom_cladelabel(node=435,label='Bt2_clade2_chimp',color='darkorange2') + 
  geom_cladelabel(node=511,label='Bt2_clade2_bonobo',color='red2') +
  geom_cladelabel(node=600,label='Bt2_clade1_chimp',color='darkorange2') + 
  geom_cladelabel(node=626,label='Bt2_clade2_gorilla',color='darkgreen') +
  #captive ape ASVs
  geom_tiplab(aes(subset=label %in% Bt2_captive_ASVs$ASV,
                  label=label),align = TRUE, color='black') +
  geom_tippoint(aes(
    subset=label %in% Bt2_captive_ASVs$ASV[Bt2_captive_ASVs$captive_gorilla_COLZ>0],
    label=label), size=1.5,x=.25,color='darkolivegreen3',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% Bt2_captive_ASVs$ASV[Bt2_captive_ASVs$captive_chimp_PC>0],
    label=label), size=1.5,x=.28,color='tan1',shape=17) 

#investigate leftover ASVs that didn't fall into in host-restricted clade
Bt2_leftover_ASVs <- Bt2_ASV_table %>% filter(cladeName=="")
table(Bt2_leftover_ASVs$HR_ASV)
table(Bt2_leftover_ASVs$ASV_sampleNum)
ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt2_tree.pdf')

#Subset tree and codivASVs to Bt3 lineage
Bt3_codiv_table <- ASV_table %>% filter(lineage == 'Bt3')
Bt3_codiv_ASVs = as.vector(Bt3_codiv_table$ASV)
Bt3_MRCA <- findMRCA(full_tree,Bt3_codiv_ASVs)
Bt3_tree <- extract.clade(full_tree,Bt3_MRCA)
Bt3_all_ASVs <- Bt3_tree$tip.label
Bt3_ASV_table <- ASV_table %>% filter(ASV %in% Bt3_all_ASVs)

Bt3_codiv_clades <- as.vector(unique(Bt3_codiv_table$codiv_clade))
(Bt3_codiv_clades <- get_codiv_table(Bt3_tree,Bt3_codiv_clades))
res <- add_HRclades_to_tree(Bt3_tree,ASV_table)
Bt3_tree_cls <- res[[1]]
Bt3_color_vec <- res[[2]]
(Bt3_captive_ASVs <- ASV_table %>% 
    filter(ASV %in% Bt3_tree$tip.label) %>% 
    filter(captive_all>0))
Bt3_captive_ASVs %>%
  select(ASV,captive_chimp_HOUZ:captive_chimp_PC) %>%
  gather(key = "host_cat", value="sample_count",-ASV) %>%
  filter(sample_count>0)

ggtree(Bt3_tree_cls, aes(color=group)) + xlim(NA, .6) +
  scale_color_manual(values=Bt3_color_vec) +
 # geom_nodepoint(aes(subset = label > .70),size=.75) + 
  geom_treescale(x=0,y=0) +
#bootstrap
  geom_cladelabel(node=625,label='Bt3_clade1_human',color='dodgerblue3',offset=.25) +
  geom_cladelabel(node=1057,label='Bt3_clade1_chimp',color='darkorange2',offset=.20) +
  geom_cladelabel(node=1077,label='Bt3_clade1_bonobo',color='red2',offset=.20) +
  geom_tiplab(aes(subset=label %in% Bt3_captive_ASVs$ASV,
                  label=label),align = TRUE, color='black') +
  geom_tippoint(aes(
    subset=label == 'ASV_934',
    label=label), size=3,color='tan1',x=.42,shape=17)  + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_bonobo_COLZ>0],
    label=label), size=3,x=.32,color='indianred2') + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_gorilla_COLZ>0],
    label=label), size=3,x=.345,color='darkolivegreen3') + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_orangutan_COLZ>0],
    label=label), size=3,x=.37,color='plum3') +
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_chimp_HOUZ>0],
    label=label), size=3,x=.395,color='tan1',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_gorilla_HOUZ>0],
    label=label), size=3,x=.42,color='darkolivegreen3',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_orangutan_HOUZ>0],
    label=label), size=3,x=.455,color='plum3',shape=15) 
#investigate leftover ASVs that didn't fall into in host-restricted clade
Bt3_leftover_ASVs <- Bt3_ASV_table %>% filter(cladeName=="")
table(Bt3_leftover_ASVs$HR_ASV)
table(Bt3_leftover_ASVs$ASV_sampleNum)

ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt3_tree.pdf')

#inset trees from Moeller et al 2016
Bt2_inset_tree <- read.tree('results/gyrb/inputs/moeller_codiv_lin_Bt2.tree')
Bt3_inset_tree <- read.tree('results/gyrb/inputs/moeller_codiv_lin_Bt3.tree')
inset_table <- read.table("results/gyrb/inputs/moeller_codiv_HRclades.txt",sep='\t',header=TRUE)

inset_table$codiv_taxa <- str_replace(inset_table$codiv_taxa, " ", "")
names(inset_table) <- c("ASV","cladeName","lineage","HR_clade")

inset_table_Bt2 <- inset_table %>% filter(ASV %in% Bt2_inset_tree$tip.label)
inset_table_Bt3 <- inset_table %>% filter(ASV %in% Bt3_inset_tree$tip.label)

Bt2_inset_tree <- keep.tip(Bt2_inset_tree,inset_table_Bt2$ASV)
Bt2_res <- add_HRclades_to_tree(Bt2_inset_tree,inset_table)
Bt2_inset_tree_cls <- Bt2_res[[1]]
(Bt2_inset_color_vec <- Bt2_res[[2]])
Bt2_inset_color_vec <- c("black","red2","darkorange2","darkgreen","red2","darkorange2")
ggtree(Bt2_inset_tree_cls, aes(color=group)) + xlim(NA, .2) +
  scale_color_manual(values=Bt2_inset_color_vec) +
  geom_treescale()
ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt2_inset_tree.pdf')

Bt3_inset_tree <- keep.tip(Bt3_inset_tree,inset_table_Bt3$ASV)
Bt3_res <- add_HRclades_to_tree(Bt3_inset_tree,inset_table)
Bt3_inset_tree_cls <- Bt3_res[[1]]
(Bt3_inset_color_vec <- Bt3_res[[2]])
Bt3_inset_color_vec <- c("black","red2","darkorange2","dodgerblue3")
ggtree(Bt3_inset_tree_cls, aes(color=group)) + xlim(NA, .4) +
  scale_color_manual(values=Bt3_inset_color_vec)  +
  geom_treescale()
ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt3_inset_tree.pdf')

