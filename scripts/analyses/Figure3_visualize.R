setwd('/Volumes/AHN/captive_ape_microbiome')
library(ape)
library(ggtree)
library(tidytree)
library(tidyverse)
library(phytools)
library(stringr)

run_Figure3 <- function(INDIR,OUTDIR) {
  full_tree <- read.tree(file.path(INDIR,'physeq_Bacteroidales.tree'))
  Bt2_inset_tree <- read.tree(file.path(INDIR,'moeller_codiv_lin_Bt2.tree'))
  Bt3_inset_tree <- read.tree(file.path(INDIR,'moeller_codiv_lin_Bt3.tree'))
  inset_table <- read.table(file.path(INDIR,"moeller_codiv_HRclades.txt"),sep='\t',header=TRUE)
  
  ASV_table <- read.table(file.path(OUTDIR,'intermediate_outputs/HRclades_Figure3_table.txt'),sep='\t',header=TRUE)
  
  
  extract_subtree <- function(tree,Bt_lin) {
    #Subset tree and codivASVs to Bt lineage
    lin_table <- ASV_table %>% filter(lineage == Bt_lin)
    lin_ASVs = as.vector(lin_table$ASV)
    MRCA <- findMRCA(tree,lin_ASVs)
    subtree <- extract.clade(tree,MRCA)
    return(subtree)
  }
  
  Bt1_tree <- extract_subtree(full_tree,'Bt1')
  ggtree(Bt1_tree) + geom_nodepoint(aes(subset = label > .50))
  
  add_HRclades_to_tree <- function(lineage_tree,table){
    #label HR clades in subtree
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
    tree_cls <- groupOTU(lineage_tree, cls)
    print(cladeHostSp)
    color_vec <- as.vector(recode(cladeHostSp, HR_wild_gorilla = "darkgreen", 
                                  HR_wild_chimp = "darkorange2",
                                  HR_wild_bonobo = "red2", 
                                  HR_human = "dodgerblue3",
                                  MX_2_wild_apes='magenta',
                                  MX_human_2_wild_apes ='magenta',
                                  MX_human_single_wild_ape='magenta'))
    (color_vec <- c("black",color_vec))
    return(c(tree_cls,color_vec))
  }
  
  Bt1_res <- add_HRclades_to_tree(Bt1_tree,ASV_table)
  Bt1_tree_cls <- Bt1_res[[1]]
  Bt1_color_vec <- Bt1_res[[2]]
  
  ggtree(Bt1_tree_cls, aes(color=group)) +
    scale_color_manual(values=Bt1_color_vec)
  
  #moeller clade labels
  unique(ASV_table$codiv_clade)
  
  get_codiv_table <- function(tree,clade){
    clade_ASVs <- ASV_table %>% filter(codiv_clade == clade)
    clade_node<- findMRCA(tree,as.vector(clade_ASVs$ASV))
    return(clade_node)
  }
  
  Bt1_gorilla_node <- get_codiv_table(Bt1_tree,"Bt1_clade1_gorilla")
  Bt1_chimp_node <- get_codiv_table(Bt1_tree,"Bt1_clade1_chimp")
  Bt1_bonobo_node <- get_codiv_table(Bt1_tree,"Bt1_clade1_bonobo")
  
  ggtree(Bt1_tree_cls, aes(color=group)) +
    scale_color_manual(values=Bt1_color_vec) +
    geom_nodepoint(aes(subset = label > .50)) +
    geom_cladelabel(node=Bt1_gorilla_node,label='Bt1_gorilla',color='green4') + 
    geom_cladelabel(node=Bt1_chimp_node,label='Bt1_chimp',color='orange') + 
    geom_cladelabel(node=Bt1_bonobo_node,label='Bt1_bonobo',color='red') 
  
  #no captive ape ASVs in Bt1 lineage
  ASV_table %>% 
    filter(ASV %in% Bt1_tree$tip.label) %>%
    filter(captive_all>0)
  
  #Subset tree and codivASVs to Bt2 lineage
  Bt2_tree <- extract_subtree(full_tree,'Bt2')
  Bt2_res <- add_HRclades_to_tree(Bt2_tree,ASV_table)
  Bt2_tree_cls <- Bt2_res[[1]]
  Bt2_color_vec <- Bt2_res[[2]]
  Bt2_clade1_gorilla_node <- get_codiv_table(Bt2_tree,"Bt2_clade1_gorilla")
  Bt2_clade1_chimp_node <- get_codiv_table(Bt2_tree,"Bt2_clade1_chimp")
  Bt2_clade1_bonobo_node <- get_codiv_table(Bt2_tree,"Bt2_clade1_bonobo")
  Bt2_clade2_chimp_node <- get_codiv_table(Bt2_tree,"Bt2_clade2_chimp")
  Bt2_clade2_bonobo_node <- get_codiv_table(Bt2_tree,"Bt2_clade2_bonobo")
  #Bt2 ASVs in captive apes 
  ASV_table %>% 
    filter(ASV %in% Bt2_tree$tip.label) %>%
    select(ASV,captive_chimp_HOUZ:captive_chimp_PC) %>%
    gather(key = "host_cat", value="sample_count",-ASV) %>%
    filter(sample_count>0)
  #Bt2 ASVs by HR_type 
  ASV_table %>% 
    filter(ASV %in% Bt2_tree$tip.label) %>%
    group_by(HR_type) %>%
    tally()
  ASV_table %>% 
    filter(ASV %in% Bt2_tree$tip.label) %>%
    group_by(ASV_HR_type) %>%
    tally()
  
  (Bt2_tree_fig <- ggtree(Bt2_tree_cls, aes(color=group)) + xlim(NA, .5) +
      scale_color_manual(values=Bt2_color_vec) + #HR clades
      geom_treescale(x=0,y=0) +
      geom_nodepoint(aes(subset = label > .50),size=.75) + #bootstrap
      #moeller codiv clades
      geom_cladelabel(node=Bt2_clade1_gorilla_node,label='Bt2_clade1_gorilla',color='darkgreen',offset=.01) + 
      geom_cladelabel(node=Bt2_clade1_chimp_node,label='Bt2_clade1_chimp',color='darkorange2',offset=.02) + 
      geom_cladelabel(node=Bt2_clade1_bonobo_node,label='Bt2_clade1_bonobo',color='red2',offset=.01) +
      geom_cladelabel(node=Bt2_clade2_chimp_node,label='Bt2_clade2_chimp',color='darkorange2',offset=.01) + 
      geom_cladelabel(node=Bt2_clade2_bonobo_node,label='Bt2_clade2_bonobo',color='red2',offset=.01) +
      #captive ape ASVs
      geom_tiplab(aes(subset=label %in% ASV_table$ASV[ASV_table$captive_all>0],
                      label=label),align = TRUE, color='black') +
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_gorilla_COLZ>0],
        label=label), size=1.5,x=.25,color='darkolivegreen3',shape=15) + 
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_chimp_PC>0],
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
        subset=label %in% ASV_table$ASV[str_detect(ASV_table$ASV_HR_type,'MX')],
        label=label), size=1,color='magenta',shape=15) + 
      theme(legend.position = "none"))
  
  ggsave(Bt2_tree_fig,filename = file.path(OUTDIR,'figures/Figure3_Bt2_tree.pdf'),width=7)
  
  
  #Subset tree and codivASVs to Bt3 lineage
  Bt3_tree <- extract_subtree(full_tree,'Bt3')
  Bt3_res <- add_HRclades_to_tree(Bt3_tree,ASV_table)
  Bt3_tree_cls <- Bt3_res[[1]]
  Bt3_color_vec <- Bt3_res[[2]]
  Bt3_clade1_human_node <- get_codiv_table(Bt3_tree,"Bt3_clade1_human")
  Bt3_clade1_chimp_node <- get_codiv_table(Bt3_tree,"Bt3_clade1_chimp")
  Bt3_clade1_bonobo_node <- get_codiv_table(Bt3_tree,"Bt3_clade1_bonobo")
  
  #Bt3 ASVs in captive apes 
  ASV_table %>% 
    filter(ASV %in% Bt3_tree$tip.label) %>%
    select(ASV,captive_chimp_HOUZ:captive_chimp_PC) %>%
    gather(key = "host_cat", value="sample_count",-ASV) %>%
    filter(sample_count>0) %>% 
    group_by(host_cat) %>% 
    tally()
  #Bt3 ASVs by HR_type 
  ASV_table %>% 
    filter(ASV %in% Bt3_tree$tip.label) %>%
    group_by(HR_type) %>%
    tally()
  ASV_table %>% 
    filter(ASV %in% Bt3_tree$tip.label) %>%
    group_by(ASV_HR_type) %>%
    tally()
  
  (Bt3_tree_fig <- ggtree(Bt3_tree_cls, aes(color=group)) + xlim(NA, .8) +
      scale_color_manual(values=Bt3_color_vec) +
      geom_nodepoint(aes(subset = label > .50),size=.75) + 
      geom_treescale(x=0,y=0) +
      geom_cladelabel(node=Bt3_clade1_human_node,label='Bt3_clade1_human',color='dodgerblue3',offset=.01) +
      geom_cladelabel(node=Bt3_clade1_chimp_node,label='Bt3_clade1_chimp',color='darkorange2',offset=.01) +
      geom_cladelabel(node=Bt3_clade1_bonobo_node,label='Bt3_clade1_bonobo',color='red2',offset=.01) +
      geom_tiplab(aes(subset=label %in%  ASV_table$ASV[ASV_table$captive_all>0],
                      label=label),align = TRUE, color='black',offset = .05) +
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_chimp_PC>0],
        label=label), size=3,color='tan1',x=.45,shape=17)  + 
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_bonobo_COLZ>0],
        label=label), size=3,x=.345,color='indianred2') + 
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_gorilla_COLZ>0],
        label=label), size=3,x=.37,color='darkolivegreen3') + 
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_orangutan_COLZ>0],
        label=label), size=3,x=.395,color='plum3') +
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_chimp_HOUZ>0],
        label=label), size=3,x=.42,color='tan1',shape=15) + 
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_gorilla_HOUZ>0],
        label=label), size=3,x=.47,color='darkolivegreen3',shape=15) + 
      geom_tippoint(aes(
        subset=label %in% ASV_table$ASV[ASV_table$captive_orangutan_HOUZ>0],
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
        subset=label %in% ASV_table$ASV[str_detect(ASV_table$ASV_HR_type,'MX')],
        label=label), size=1,color='magenta',shape=15) + 
      theme(legend.position = "none"))
  
  ggsave(Bt3_tree_fig,filename = file.path(OUTDIR,'figures/FigureSup_Bt3_tree.pdf'),width=7)
  
}

INDIR = 'results/gyrb/inputs/ps_Bacteroidales_asvRefTree'
OUTDIR = 'results/gyrb/analyses/ps_Bacteroidales_asvRefTree'
run_Figure3(INDIR,OUTDIR)
INDIR = 'results/gyrb/inputs/ps_Bacteroidales_asvTree'
OUTDIR = 'results/gyrb/analyses/ps_Bacteroidales_asvTree'
run_Figure3(INDIR,OUTDIR)
