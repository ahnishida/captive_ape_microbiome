library(ape)
library(ggtree)
library(tidyverse)
library(phytools)
library(cowplot)
library(ggplot2)

setwd('/Volumes/AHN/captive_ape_microbiome') #SET WORKING DIR

#summarize 16S and gyrB ASV distrbutions
HR_16S <- read.table('results/16s/analyses/tables/16S_ASVs_summary.txt',sep='\t',header=T)
HR_16S <- HR_16S %>% filter(Order == 'Bacteroidales') 
All_16S <- HR_16S %>% group_by(HR_type) %>% tally()
All_16S$cat <- '16S Bacteroidales ASVs all samples'
CP_16S <- HR_16S %>% filter(CP_pres=='True') %>% group_by(HR_type) %>% tally()
CP_16S$cat <- '16S Bacteroidales ASVs captive apes'

HR_gyrb <- read.table('results/gyrb/analyses/tables/gyrb_asv_hr_table.txt',sep='\t',header=T)
All_gyrb <- HR_gyrb  %>% group_by(HR_type) %>% tally()
All_gyrb$cat <- 'gyrb ASVs all samples'
CP_gyrb <- HR_gyrb %>% filter(CP_pres=='True') %>% group_by(HR_type) %>% tally()
CP_gyrb$cat <- 'gyrb ASVs captive apes'
HR_16S_gyrb <- bind_rows(All_16S,CP_16S,All_gyrb,CP_gyrb)
print('summary of 16S and gyrb ASVs observed')
HR_16S_gyrb %>% group_by(cat) %>% summarise(sum(n))

HRpalette <- as.vector(recode(unique(HR_16S_gyrb$HR_type), HR_wild_gorilla = "darkgreen", 
                               HR_wild_chimp = "darkorange2",
                               HR_wild_bonobo = "red2", 
                               HR_human = "dodgerblue",
                               MX_human_wild_apes = "maroon",
                               MX_wild_apes = "goldenrod", 
                               Unique_CP='purple'))
plot_HR_16S_gyrb <- ggplot(HR_16S_gyrb, aes(fill=HR_type, y=n, x=cat)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=HRpalette) + 
  theme_cowplot() +  
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  theme(axis.title.x=element_blank()) 
ggsave(plot_HR_16S_gyrb,file='results/gyrb/analyses/figures/FigureS9_HR_16S_gyrb.pdf')

#FIGURE 4
tree_file <-  file.path('results/gyrb/analyses/intermediate_outputs/HRclades_wholetree.tre')
clades_table <- file.path('results/gyrb/analyses/intermediate_outputs/HRclades_wholetree_table.txt')

#read in clades data
clades <- data.frame(read.table(clades_table,sep='\t',header=TRUE))

print('Summary of host restricted gyrB clades')
print('gyrB ASVs distributed among host-restricted, mixed host, and unique to captive clades')
(clades %>% group_by(HR_cat) %>% summarise(sum(ASVsNum)))
print('clades across all samples')
(c1_summary <- clades %>% filter(heatmap_col1 !='Blank') %>% group_by(heatmap_col1) %>% tally())
print(c(sum(c1_summary$n),'total clades'))
print('clades present in 25% of individuals of any host species in captivity or wild')
(c2_summary <- clades %>% filter(heatmap_col2 !='Blank') %>% group_by(heatmap_col2) %>% tally())
print(c(sum(c2_summary$n),'total clades'))
print('clades present in 25% of individuals of any host species in captivity')
(c3_summary <- clades %>% filter(heatmap_col3 !='Blank') %>% group_by(heatmap_col3) %>% tally())
print(c(sum(c3_summary$n),'total clades'))

#subset to clades that are observed in threshold cutoff of the individuals with one captive ape species
threshold <- .25
captive_clades <- clades %>%
    filter(heatmap_col3 !='Blank')

#reformat 
captive_clades_long <- captive_clades %>% 
  select(cladeName,wild_gorilla,wild_chimp,wild_bonobo,industrialized_human,non_industrialized_human,captive_chimp,captive_gorilla,captive_orangutan,captive_bonobo) %>%
  gather(key = "host_cat", value="percent_samples",-cladeName)

#reorder sample description 
new_order = c("non_industrialized_human","industrialized_human","wild_bonobo", "wild_chimp","wild_gorilla",
              "captive_bonobo","captive_chimp","captive_gorilla","captive_orangutan")
captive_clades_long$host_cat <- factor(captive_clades_long$host_cat, levels = new_order)
captive_clades_long <- captive_clades_long %>% filter(percent_samples>0)
#read in phylogeny
full_tree <- ape::read.tree(tree_file)

#add moeller clades into the tree
Bt1_lineage_clades = clades %>% filter(lineage == 'Bt1')
Bt1_lineage_clades$cladeName
Bt1_node = findMRCA(full_tree, as.vector(Bt1_lineage_clades$cladeName))
Bt2_lineage_clades = clades %>% filter(lineage == 'Bt2')
Bt2_lineage_clades$cladeName
Bt2_node = findMRCA(full_tree, as.vector(Bt2_lineage_clades$cladeName))
Bt3_lineage_clades = clades %>% filter(lineage == 'Bt3')
Bt3_lineage_clades = intersect(Bt3_lineage_clades$cladeName,full_tree$tip.label)
Bt3_node = findMRCA(full_tree, as.vector(Bt3_lineage_clades))

#tree plot
tree_plot <- ggtree(full_tree,ladderize = F)  %<+% clades + 
  geom_cladelabel(node=Bt1_node,label='Bt1',color='magenta') + 
  geom_cladelabel(node=Bt2_node,label='Bt2',color='magenta') +
  geom_cladelabel(node=Bt3_node,label='Bt3',color='magenta') +
  geom_treescale()

#add heatmap
clades_heatmap <- clades %>% 
  column_to_rownames("cladeName") %>% 
  select(heatmap_col1,heatmap_col2,heatmap_col3)

(tree_heatmap <-gheatmap(tree_plot, clades_heatmap, offset = .01, width=0.5) + 
  scale_fill_manual(values=c("white","dodgerblue3","red2","darkorange2","darkgreen","chocolate4","purple")))
ggsave(tree_heatmap,filename = 'results/gyrb/analyses/figures/Figure4_tree_plot.pdf',width=5,height=10)


#Dotplots

#generate color scales
dotplot_color_scale <- levels(droplevels(factor(captive_clades_long$host_cat)))
(dotplot_color_scale <- as.vector(recode(dotplot_color_scale, wild_gorilla = "darkgreen", 
                                         wild_chimp = "darkorange2",
                                         wild_bonobo = "red2", 
                                         industrialized_human = "blue",
                                         non_industrialized_human="skyblue2",
                                         captive_gorilla = "darkolivegreen3", 
                                         captive_chimp = "tan1",
                                         captive_bonobo = "indianred2", 
                                         captive_orangutan="plum3")))

HRclade_color_scale <- levels(droplevels(factor(captive_clades$HR_type)))
(HRclade_color_scale <- as.vector(recode(HRclade_color_scale, HR_wild_gorilla = "darkgreen", 
                                         HR_wild_chimp = "darkorange2",
                                         HR_wild_bonobo = "red2", 
                                         HR_human = "dodgerblue3",
                                         MX_human_wild_apes = "chocolate4",
                                         MX_wild_apes = "chocolate4", 
                                         Unique_CP = "purple")))

#change order of taxa in dotplot
ordered_tips = full_tree$tip.label[full_tree$tip.label %in% as.character(captive_clades$cladeName)]
setdiff(ordered_tips,as.character(unique(captive_clades_long$cladeName)))
setdiff(as.character(unique(captive_clades_long$cladeName)),ordered_tips)
captive_clades_long$cladeName = factor(captive_clades_long$cladeName,levels = ordered_tips)
captive_clades$cladeName = factor(captive_clades$cladeName,levels = ordered_tips)

#generate dot plot
dotplot <- ggplot(captive_clades_long)  +
  geom_point(aes(x=host_cat,y=cladeName,color = host_cat,size=percent_samples))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank()) +
  scale_color_manual(values=dotplot_color_scale)
dotplot

#taxonomy_dotplot 
(taxonomy_dotplot <- ggplot(captive_clades) +
  geom_point(aes(x=0,y=cladeName,color = cladeTax))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank()))

#add HR_type
(HR_type_dotplot <- ggplot(captive_clades) +
  geom_point(aes(x=0,y=cladeName,color = HR_type))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank()) +
  scale_color_manual(values=HRclade_color_scale))

#put all the dotplots together
all_dotplots <- plot_grid(HR_type_dotplot,taxonomy_dotplot,dotplot,nrow = 1)
ggsave(all_dotplots,filename = 'results/gyrb/analyses/figures/Figure4_all_dotplots.pdf',height=10,width=18)

#SUBTREE FIGURE 5, S10 - sub

full_tree <- read.tree(file.path('results/gyrb/inputs/physeq_Bacteroidales_ASVs_ref.tree'))
Bt2_inset_tree <- read.tree(file.path('results/gyrb/inputs/moeller_codiv_lin_Bt2.tree'))
Bt3_inset_tree <- read.tree(file.path('results/gyrb/inputs/moeller_codiv_lin_Bt3.tree'))
inset_table <- read.table(file.path('results/gyrb/inputs/moeller_codiv_HRclades.txt'),sep='\t',header=TRUE)
ASV_table <- read.table(file.path('results/gyrb/analyses/intermediate_outputs/HRclades_subtrees_table.txt'),sep='\t',header=TRUE)

print('distribution of ASVs in clades')

  
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
      #print(clade)
      (clade_table <- table %>% filter(cladeName == clade))
      (clade_ASVs <- as.character(clade_table$ASV))
      cls[[length(cls)+1]] <- clade_ASVs
    }
    tree_cls <- groupOTU(lineage_tree, cls)
    #print(cladeHostSp)
    color_vec <- as.vector(recode(cladeHostSp, HR_wild_gorilla = "darkgreen", 
                                  HR_wild_chimp = "darkorange2",
                                  HR_wild_bonobo = "red2", 
                                  HR_human = "dodgerblue3",
                                  MX_wild_apes='magenta',
                                  MX_human_wild_apes ='magenta'))
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
  
ggsave(Bt2_tree_fig,filename = file.path('results/gyrb/analyses/figures/Figure5_Bt2_tree.pdf'),width=5)
  
  
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
  
ggsave(Bt3_tree_fig,filename = file.path('results/gyrb/analyses/figures/FigureS10_Bt3_tree.pdf'),width=7)
  

