setwd('/Volumes/AHN/captive_ape_microbiome')
library(ape)
library(ggtree)
library(tidytree)
library(tidyverse)
library(ggrepel)
library(phytools)
library(gggibbous)
library(aplot)

#inputs/outputs
tree_file <- 'results/gyrb/analyses/codiv_moeller_ASVs/full_tree_collapsed.tre'
clades_table <- 'results/gyrb/analyses/codiv_moeller_ASVs/full_tree_clades_collapsed_table.txt'
plot_file <- 'results/gyrb/analyses/codiv_moeller_ASVs/Figure3.pdf'

#read in clades data
clades <- data.frame(read.table(clades_table,sep='\t',header=TRUE))
#subset to clades that are observed in 25% of the individuals with one captive ape species
(captive_clades <- clades %>%
    filter(captive_chimp>.25|captive_gorilla>.25|captive_orangutan>.25|captive_bonobo>.25))
#reformat so gggibbous can plot
captive_clades_long <- captive_clades %>% 
  select(cladeName,wild_gorilla,wild_chimp,wild_bonobo,western_human,non_western_human,captive_chimp,captive_gorilla,captive_orangutan,captive_bonobo) %>%
  gather(key = "host_cat", value="percent_samples",-cladeName)
#create color scale for moon plot
unique(captive_clades_long$host_cat)
moon_color_scale <- levels(as.factor(captive_clades_long$host_cat))
(moon_color_scale <- as.vector(recode(moon_color_scale, wild_gorilla = "darkgreen", 
                                 wild_chimp = "darkorange2",
                                 wild_bonobo = "red2", 
                                 western_human = "blue",
                                 non_western_human="skyblue2",
                                 captive_gorilla = "darkolivegreen3", 
                                 captive_chimp = "tan1",
                                 captive_bonobo = "indianred2", 
                                 captive_orangutan="plum3")))
#generate moon plot
(moon_plot <- ggplot(captive_clades_long)  +
    geom_moon(aes(ratio = percent_samples, x=host_cat,  y=cladeName, color = host_cat, fill = host_cat),right =FALSE) +
    geom_moon(aes(ratio = 1 - percent_samples, x=host_cat,  y=cladeName, color = host_cat), right = TRUE) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank()) +
    scale_color_manual(values=moon_color_scale) +
    scale_fill_manual(values=moon_color_scale)) 

#read in phylogeny
full_tree <- read.tree(tree_file)

#determine clades present in 1% of individuals of any sample group
perc_cutoff=.01
(clades_over_perc_cutoff <- clades %>%
    filter(captive_chimp>perc_cutoff|captive_gorilla>perc_cutoff|captive_orangutan>perc_cutoff|captive_bonobo>perc_cutoff|
           wild_chimp>perc_cutoff|wild_bonobo>perc_cutoff|wild_gorilla>perc_cutoff|
           western_human>perc_cutoff|non_western_human>perc_cutoff))

#code in tree color scale
(tree_color_scale <- levels(as.factor(clades$HR_clade)))

(tree_color_scale <- as.vector(recode(tree_color_scale, wild_gorilla = "darkgreen", 
                                wild_chimp = "darkorange2",
                                wild_bonobo = "red2", 
                                human = "dodgerblue3",
                                mixed="chocolate4",
                                captive="purple")))

(tree_plot <- ggtree(full_tree,ladderize = F) %<+% clades + 
  geom_tippoint(aes(color=HR_clade,subset=label %in% as.character(clades_over_perc_cutoff$cladeName)),size=1) + 
  scale_color_manual(values=tree_color_scale) +
  geom_tiplab(aes(subset=label %in% as.character(captive_clades$cladeName)),
                  align=TRUE) + xlim(NA, 1)+ 
  theme(legend.position='bottom'))

moon_plot %>% insert_left(tree_plot) 

#change order of taxa in moon plot
ordered_tips = full_tree$tip.label[full_tree$tip.label %in% as.character(captive_clades$cladeName)]
captive_clades_long$cladeName = factor(captive_clades_long$cladeName,levels = ordered_tips)
(moon_plot <- ggplot(captive_clades_long)  +
    geom_moon(aes(ratio = percent_samples, x=host_cat,  y=cladeName, color = host_cat, fill = host_cat),right =FALSE) +
    geom_moon(aes(ratio = 1 - percent_samples, x=host_cat,  y=cladeName, color = host_cat), right = TRUE) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank()) +
    scale_color_manual(values=moon_color_scale) +
    scale_fill_manual(values=moon_color_scale)) 
moon_plot %>% insert_left(tree_plot) 
head(clades)
#add taxonomy to tree
ggtree(full_tree,ladderize = F) %<+% clades + 
  geom_tippoint(aes(color=cladeTax)) + 
  theme(legend.position='bottom')
#define nodes for genera
table(captive_clades$Gen)
Bacteroides_clades <- clades %>% filter(cladeTax =='Bacteroidaceae_Bacteroides')
Bacteroides_clades <- Bacteroides_clades$cladeName
(Bacteroides_clades <- as.character(Bacteroides_clades))
Bacteroides_node  <- findMRCA(full_tree,Bacteroides_clades)
Prevotella_clades <- clades %>% filter(cladeTax =='Bacteroidaceae_Prevotella')
Prevotella_clades <- Prevotella_clades$cladeName
Prevotella_clades <- as.character(Prevotella_clades)
Prevotella_node  <- findMRCA(full_tree,Prevotella_clades)
Tannerellaceae_clades <- clades %>% filter(cladeTax =='Tannerellaceae_Parabacteroides')
Tannerellaceae_clades <- Tannerellaceae_clades$cladeName
Tannerellaceae_clades <- as.character(Tannerellaceae_clades)
Tannerellaceae_node  <- findMRCA(full_tree,Tannerellaceae_clades)

(tree_plot <- ggtree(full_tree,ladderize = F) %<+% clades + 
      geom_tippoint(aes(color=HR_clade,subset=label %in% as.character(clades_over_perc_cutoff$cladeName)),size=1) + 
      geom_hilight(node=Bacteroides_node, fill="lightgreen",alpha=.25) +
      geom_hilight(node=Prevotella_node, fill="gold",alpha=.25) +
      geom_hilight(node=Tannerellaceae_node , fill="lightblue",alpha=.25) +
 
     scale_color_manual(values=tree_color_scale) +
     geom_tiplab(aes(color= HR_clade,
                     subset=label %in% as.character(captive_clades$cladeName)),
                 align=TRUE) + xlim(NA, 1)+ 
     theme(legend.position='bottom')
    )
(tree_moon_plot <- moon_plot %>% insert_left(tree_plot))
ggsave(tree_moon_plot,width=15,height=10,file=plot_file)
