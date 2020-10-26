print(.libPaths()) 
setwd('/Volumes/AHN/captive_ape_microbiome')
library(ape)
library(ggtree)
library(tidyverse)
library(phytools)
library(cowplot)
library(ggplot2)

run_Figure2 <- function(INDIR,threshold){
#inputs
tree_file <-  file.path(INDIR,'intermediate_outputs/HRclades_Figure2.tre')
clades_table <- file.path(INDIR,'intermediate_outputs/HRclades_Figure2_table.txt')

#read in clades data
clades <- data.frame(read.table(clades_table,sep='\t',header=TRUE))
#subset to clades that are observed in threshold cutoff of the individuals with one captive ape species
captive_clades <- clades %>%
    filter(captive_chimp>threshold|captive_gorilla>threshold|captive_orangutan>threshold|captive_bonobo>threshold)
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
ggsave(tree_heatmap,filename=file.path(INDIR,'figures/Figure2_tree_plot.pdf'),width=5,height=10)

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
                                         MX_2_wild_apes = "chocolate4",
                                         MX_human_2_wild_apes = "chocolate4", 
                                         MX_human_single_wild_ape = "chocolate4",
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
ggsave(all_dotplots,filename=file.path(INDIR,'figures/Figure2_all_dotplots.pdf'),height=10,width=18)

}
run_Figure2('results/gyrb/analyses/ps_Bacteroidales_asvTree',.15)
run_Figure2('results/gyrb/analyses/ps_Bacteroidales_asvRefTree',.15)
