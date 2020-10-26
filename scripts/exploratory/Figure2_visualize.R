library(ape)
library(ggtree)
library(tidyverse)
library(phytools)
library(cowplot)
library(ggplot2)

setwd('/Volumes/AHN/captive_ape_microbiome/results/gyrb/exploratory/gyrb_amp_meta_datasets_10.2')

run_Figure2 <- function(folder){

#inputs
tree_file <-  file.path(folder,'HRclades_Figure2.tre')
clades_table <- file.path(folder,'HRclades_Figure2_table.txt')

#read in clades data
clades <- data.frame(read.table(clades_table,sep='\t',header=TRUE))
#subset to clades that are observed in 25% of the individuals with one captive ape species
threshold <- .15
(captive_clades <- clades %>%
    filter(captive_chimp>threshold|captive_gorilla>threshold|captive_orangutan>threshold|captive_bonobo>threshold))
captive_clades
#reformat 
captive_clades_long <- captive_clades %>% 
  select(cladeName,wild_gorilla,wild_chimp,wild_bonobo,western_human,non_western_human,captive_chimp,captive_gorilla,captive_orangutan,captive_bonobo) %>%
  gather(key = "host_cat", value="percent_samples",-cladeName)
#create color scale for moon plot
unique(captive_clades_long$host_cat)
new_order = c("non_western_human","western_human","wild_bonobo", "wild_chimp","wild_gorilla",
                                                     "captive_bonobo","captive_chimp","captive_gorilla","captive_orangutan")
captive_clades_long$host_cat <- factor(captive_clades_long$host_cat, levels = new_order)
captive_clades_long <- captive_clades_long %>% filter(percent_samples>0)


dotplot_color_scale <- levels(as.factor(captive_clades_long$host_cat))
(dotplot_color_scale <- as.vector(recode(dotplot_color_scale, wild_gorilla = "darkgreen", 
                                 wild_chimp = "darkorange2",
                                 wild_bonobo = "red2", 
                                 western_human = "blue",
                                 non_western_human="skyblue2",
                                 captive_gorilla = "darkolivegreen3", 
                                 captive_chimp = "tan1",
                                 captive_bonobo = "indianred2", 
                                 captive_orangutan="plum3")))

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

#read in phylogeny
full_tree <- read.tree(tree_file)

#code in tree color scale
(tree_color_scale <- levels(as.factor(clades$HR_type)))
(tree_color_scale <- as.vector(recode(tree_color_scale, HR_wild_gorilla = "darkgreen", 
                                HR_wild_chimp = "darkorange2",
                                HR_wild_bonobo = "red2", 
                                HR_human = "dodgerblue3",
                                MX_2_wild_apes="chocolate4",
                                MX_3_wild_apes="chocolate4",
                                MX_human_single_wild_ape="chocolate4",
                                MX_human_2_wild_apes="chocolate4",
                                Unique_CP="purple")))

tree_plot <- ggtree(full_tree,ladderize = F)  %<+% clades + 
  geom_tippoint(aes(color=HR_type),size=1) +
  geom_tiplab(aes(subset=label %in% as.character(captive_clades$cladeName)), align=TRUE) + 
  xlim(NA, 1)+ 
  scale_color_manual(values=tree_color_scale) +
  theme(legend.position='bottom')
tree_plot <- ggtree(full_tree,ladderize = F)  %<+% clades + 
  geom_tippoint(aes(color=HR_type),size=1) +
  geom_tiplab(aes(subset=label %in% as.character(captive_clades$cladeName)), align=TRUE) +
  scale_color_manual(values=tree_color_scale) +
  theme(legend.position='bottom')+
  xlim(NA, 1.5)
tree_plot

#heatmap
tree_plot <- ggtree(full_tree,ladderize = F) 
clades_heatmap <- clades %>% 
  column_to_rownames("cladeName") %>% 
  select(heatmap_col1,heatmap_col2,heatmap_col3)
clades_heatmap 
gheatmap(tree_plot, clades_heatmap, offset = .01, width=0.5) + 
  scale_fill_manual(values=c("white","dodgerblue3","red2","darkorange2","darkgreen","chocolate4","purple"))

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


tree_plot <- ggtree(full_tree,ladderize = F)  %<+% clades + 
  geom_cladelabel(node=Bt1_node,label='Bt1',color='magenta') + 
  geom_cladelabel(node=Bt2_node,label='Bt2',color='magenta') +
  geom_cladelabel(node=Bt3_node,label='Bt3',color='magenta') +
  geom_treescale()
tree_plot 
tree_heatmap <- gheatmap(tree_plot, clades_heatmap, offset = .015, width=0.5) + 
  scale_fill_manual(values=c("white","dodgerblue3","red2","darkorange2","darkgreen","chocolate4","purple"))

ggsave(tree_heatmap,filename=file.path(folder,'tree_plot.pdf'),width=5,height=10)


#see what they look like together
plot_grid(tree_plot,dotplot)

#change order of taxa in dotplot
ordered_tips = full_tree$tip.label[full_tree$tip.label %in% as.character(captive_clades$cladeName)]
captive_clades_long$cladeName = factor(captive_clades_long$cladeName,levels = ordered_tips)
(dotplot <- ggplot(captive_clades_long)  +
    geom_point(aes(x=host_cat,y=cladeName,color = host_cat,size=percent_samples))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank()) +
    scale_color_manual(values=dotplot_color_scale))
dotplot_color_scale_sansbonbo = c("skyblue2","blue","darkorange2","darkgreen",
  "indianred2","tan1","darkolivegreen3","plum3")
#taxonomy_dotplot 
captive_clades$cladeName = factor(captive_clades$cladeName,levels = ordered_tips)
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
  scale_color_manual(values=c("dodgerblue3","darkorange2","chocolate4",
                              "chocolate4","chocolate4","purple"))) 

#put all the dotplots together
all_dotplots <- plot_grid(HR_type_dotplot,taxonomy_dotplot,dotplot,nrow = 1)

ggsave(file.path(folder,'all_dotplots.pdf'),height=10,width=18)

}

run_Figure2('BtRef_AlignMafft')
run_Figure2('BtRef_AlignTrans')
run_Figure2('GemRef_AlignMafft')
run_Figure2('GemRef_AlignTrans')
run_Figure2('noRef_AlignMafft')
run_Figure2('Bt_Family5_ref')
run_Figure2('noRef_AlignTrans')
run_Figure2('BacteroidalesOnly_AlignTrans')

