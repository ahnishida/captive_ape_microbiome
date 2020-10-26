
library(ape)
library(ggtree)
library(tidytree)
library(tidyverse)
library(phytools)
setwd('/Volumes/AHN/captive_ape_microbiome/results/gyrb/processing/gyrb_amp_meta_datasets_10.2/')

folder <- 'BtRef_AlignMafft'
run_Figure3 <- function(folder) {

full_tree <- read.tree(file.path(folder,'ASVs_filtered_ref_aln_trim_rooted.tre'))
ASV_table <- read.table(file.path(folder,"Figure3_table.txt"),sep='\t',header=TRUE)
Bt2_inset_tree <- read.tree('moeller_codiv_lin_Bt2.tree')
Bt3_inset_tree <- read.tree('moeller_codiv_lin_Bt3.tree')
inset_table <- read.table("moeller_codiv_HRclades.txt",sep='\t',header=TRUE)

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
  geom_cladelabel(node=76,label='Bt1_gorilla',color='green4') + 
  geom_cladelabel(node=46,label='Bt1_chimp',color='orange') + 
  geom_cladelabel(node=53,label='Bt1_bonobo',color='red') 

#import in HRclades filter to those present in lineage
add_HRclades_to_tree <- function(lineage_tree,table){
  (taxa<-lineage_tree$tip.label) 
  (lineage_table <- table %>% 
      dplyr::filter(ASV %in% taxa & cladeName != "") %>% 
      dplyr::group_by(cladeName) %>% 
      dplyr::sample_n(1))
  (cladeNames = as.vector(lineage_table$cladeName))
  (cladeHostSp = as.vector(lineage_table$HR_type))
  cls <- list()
  for (clade in cladeNames) {
    print(clade)
    (clade_table <- table %>% filter(cladeName == clade))
    (clade_ASVs <- as.character(clade_table$ASV))
    cls[[length(cls)+1]] <- clade_ASVs
  }
  cls
  tree_cls <- groupOTU(lineage_tree, cls)
  print(cladeHostSp)
  color_vec <- as.vector(recode(cladeHostSp, HR_wild_gorilla = "darkgreen", 
                                HR_wild_chimp = "darkorange2",
                                HR_wild_bonobo = "red2", 
                                HR_human = "dodgerblue3",
                                MX_2_wild_apes='magenta',
                                MX_human_2_wild_apes ='darkmagenta',
                                MX_human_single_wild_ape='darkmagenta'))
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
  geom_cladelabel(node=76,label='Bt1 gorilla',color='darkgreen') + 
  geom_cladelabel(node=46,label='Bt1 chimp',color='darkorange2') + 
  geom_cladelabel(node=53,label='Bt1 bonobo',color='red2') +
  scale_color_manual(values=color_vec)
#ggsave('results/gyrb/analyses/codiv_moeller_ASVs/Bt1_tree.pdf')

#no captive ape ASVs in lineage
Bt1_ASV_table %>%
  filter(captive_all>0)

#practice adding geomtip colors
ggtree(tree_cls) + xlim(NA, .3) +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_wild_gorilla'],
    label=label), size=1.5,color='tan1',shape=17) 



#Subset tree and codivASVs to Bt2 lineage
Bt2_codiv_table <- ASV_table %>% filter(lineage == 'Bt2')
Bt2_codiv_ASVs = as.vector(Bt2_codiv_table$ASV)
Bt2_MRCA <- findMRCA(full_tree,Bt2_codiv_ASVs)
Bt2_tree <- extract.clade(full_tree,Bt2_MRCA)
Bt2_all_ASVs <- Bt2_tree$tip.label
Bt2_ASV_table <- ASV_table %>% filter(ASV %in% Bt2_all_ASVs)
print(unique(Bt2_ASV_table$HR_type))
print(unique(Bt2_ASV_table$ASV_HR_type))
#add moeller codiv clades 
Bt2_codiv_clades <- as.vector(unique(Bt2_codiv_table$codiv_clade))
(Bt2_codiv_clades <- get_codiv_table(Bt2_tree,Bt2_codiv_clades))

#define host-restricted and mixed clades
res <- add_HRclades_to_tree(Bt2_tree,ASV_table)
Bt2_tree_cls <- res[[1]]
Bt2_color_vec <- res[[2]]
Bt2_color_vec
#determine whether captive apes harbor any ASVs in the tree
Bt2_captive_ASVs <- Bt2_ASV_table %>%
  filter(captive_all>0)

ggtree(Bt2_tree_cls, aes(color=group)) +
  scale_color_manual(values=Bt2_color_vec)

ggtree(Bt2_tree_cls, aes(color=group)) + xlim(NA, .5) +
  scale_color_manual(values=Bt2_color_vec) + #HR clades
  geom_treescale(x=0,y=0) +
  #geom_nodepoint(aes(subset = label > .70),size=.75) + #bootstrap
  #moeller codiv clades
  #geom_cladelabel(node=502,label='Bt2_clade2_bonobo',color='red2') + 
  #geom_cladelabel(node=536,label='Bt2_clade2_chimp',color='darkorange2') + 
  #geom_cladelabel(node=581,label='Bt2_clade2_bonobo',color='red2') +
  #geom_cladelabel(node=576,label='Bt2_clade1_chimp',color='darkorange2') + 
  #geom_cladelabel(node=713,label='Bt2_clade2_gorilla',color='darkgreen') +
  #captive ape ASVs
  geom_tiplab(aes(subset=label %in% Bt2_captive_ASVs$ASV,
                  label=label),align = TRUE, color='black') +
  geom_tippoint(aes(
    subset=label %in% Bt2_captive_ASVs$ASV[Bt2_captive_ASVs$captive_gorilla_COLZ>0],
    label=label), size=1.5,x=.25,color='darkolivegreen3',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% Bt2_captive_ASVs$ASV[Bt2_captive_ASVs$captive_chimp_PC>0],
    label=label), size=1.5,x=.28,color='tan1',shape=17) +
  #HR_ASVs labels
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_wild_chimp'],
    label=label), size=1,color='darkorange2',shape=15) +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_wild_bonobo'],
    label=label), size=1,color='red2',shape=15) +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_wild_gorilla'],
    label=label), size=1,color='darkgreen',shape=15) +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_human'],
    label=label), size=1,color='dodgerblue3',shape=15)  +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='MX_2_wild_apes'],
    label=label), size=1,color='darkmagenta',shape=15)

#investigate leftover ASVs that didn't fall into in host-restricted clade
Bt2_leftover_ASVs <- Bt2_ASV_table %>% filter(cladeName=="")
table(Bt2_leftover_ASVs$HR_ASV)
table(Bt2_leftover_ASVs$ASV_sampleNum)
ggsave(file.path(folder,'Bt2_tree.pdf'),width=4)

#Subset tree and codivASVs to Bt3 lineage
Bt3_codiv_table <- ASV_table %>% filter(lineage == 'Bt3')
Bt3_codiv_ASVs = as.vector(Bt3_codiv_table$ASV)
Bt3_MRCA <- findMRCA(full_tree,Bt3_codiv_ASVs)
Bt3_tree <- extract.clade(full_tree,Bt3_MRCA)
Bt3_all_ASVs <- Bt3_tree$tip.label
Bt3_ASV_table <- ASV_table %>% filter(ASV %in% Bt3_all_ASVs)
print(unique(Bt3_ASV_table$HR_type))
print(unique(Bt3_ASV_table$ASV_HR_type))
plot(Bt3_tree)
#why are some tree tips black? they are ref seqs
length(Bt3_ASV_table$ASV)
length(Bt3_tree$tip.label)
setdiff(Bt3_tree$tip.label,Bt3_ASV_table$ASV)


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
ggtree(Bt3_tree_cls, aes(color=group)) + xlim(NA, .8) +
  scale_color_manual(values=Bt3_color_vec) 
ggtree(Bt3_tree_cls, aes(color=group)) + xlim(NA, .8) +
  scale_color_manual(values=Bt3_color_vec) +
 # geom_nodepoint(aes(subset = label > .70),size=.75) + 
  geom_treescale(x=0,y=0) +
#bootstrap
  #geom_cladelabel(node=906,label='Bt3_clade1_human',color='dodgerblue3',offset=.30) +
  #geom_cladelabel(node=1394,label='Bt3_clade1_chimp',color='darkorange2',offset=.25) +
  #geom_cladelabel(node=1406,label='Bt3_clade1_bonobo',color='red2',offset=.25) +
  geom_tiplab(aes(subset=label %in% Bt3_captive_ASVs$ASV,
                  label=label),align = TRUE, color='black',offset = .05) +
  geom_tippoint(aes(
    subset=label == 'ASV_934',
    label=label), size=3,color='tan1',x=.42,shape=17)  + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_bonobo_COLZ>0],
    label=label), size=3,x=.345,color='indianred2') + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_gorilla_COLZ>0],
    label=label), size=3,x=.37,color='darkolivegreen3') + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_orangutan_COLZ>0],
    label=label), size=3,x=.395,color='plum3') +
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_chimp_HOUZ>0],
    label=label), size=3,x=.42,color='tan1',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_gorilla_HOUZ>0],
    label=label), size=3,x=.47,color='darkolivegreen3',shape=15) + 
  geom_tippoint(aes(
    subset=label %in% Bt3_captive_ASVs$ASV[Bt3_captive_ASVs$captive_orangutan_HOUZ>0],
    label=label), size=3,x=.495,color='plum3',shape=15) +
  #HR_ASVs labels
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_wild_chimp'],
    label=label), size=1,color='darkorange2',shape=15) +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_wild_bonobo'],
    label=label), size=1,color='red2',shape=15) +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_wild_gorilla'],
    label=label), size=1,color='darkgreen',shape=15) +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='HR_human'],
    label=label), size=1,color='dodgerblue3',shape=15)  +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='MX_2_wild_apes'],
    label=label), size=1,color='magenta',shape=15)  +
  geom_tippoint(aes(
    subset=label %in% ASV_table$ASV[ASV_table$ASV_HR_type=='MX_human_single_wild_ape'],
    label=label), size=1,color='darkmagenta',shape=15)  

#investigate leftover ASVs that didn't fall into in host-restricted clade
Bt3_leftover_ASVs <- Bt3_ASV_table %>% filter(cladeName=="")
table(Bt3_leftover_ASVs$HR_ASV)
table(Bt3_leftover_ASVs$ASV_sampleNum)

ggsave(file.path(folder,'Bt3_tree.pdf'),width=7)
}