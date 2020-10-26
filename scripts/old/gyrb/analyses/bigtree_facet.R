library(ape)
library(ggtree)
library(tidytree)
library(tidyverse)
library(ggrepel)
library(phytools)
library(gggibbous)
library(aplot)

setwd('/Volumes/AHN/captive_ape_microbiome/results/gyrb/analyses/figures')

#inputs/outputs
tree_file <-  'HRclades_Figure2.tre'
clades_table <-'HRclades_Figure2_table.txt'
plot_file <- 'Figure2_raw.pdf'

#read in clades data
clades <- data.frame(read.table(clades_table,sep='\t',header=TRUE))
#subset to clades that are observed in 25% of the individuals with one captive ape species
threshold <- .15
(captive_clades <- clades %>%
    filter(captive_chimp>threshold|captive_gorilla>threshold|captive_orangutan>threshold|captive_bonobo>threshold))
#reformat 
captive_clades_long <- captive_clades %>% 
  select(cladeName,wild_gorilla,wild_chimp,wild_bonobo,western_human,non_western_human,captive_chimp,captive_gorilla,captive_orangutan,captive_bonobo) %>%
  gather(key = "host_cat", value="percent_samples",-cladeName)
#create color scale for moon plot
unique(captive_clades_long$host_cat)
captive_clades_long$host_cat <- factor(captive_clades_long$host_cat)
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
    scale_color_manual(values=c("skyblue2","blue","red2","darkorange2","darkgreen",
                               "indianred2","tan1","darkolivegreen3","plum3"))
dotplot
#read in phylogeny
full_tree <- read.tree(tree_file)

#determine clades present in 1% of individuals of any sample group
perc_cutoff=.00
(clades_over_perc_cutoff <- clades %>%
    filter(captive_chimp>perc_cutoff|captive_gorilla>perc_cutoff|captive_orangutan>perc_cutoff|captive_bonobo>perc_cutoff|
           wild_chimp>perc_cutoff|wild_bonobo>perc_cutoff|wild_gorilla>perc_cutoff|
           western_human>perc_cutoff|non_western_human>perc_cutoff))

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
tree_plot               

library(cowplot)
plot_grid(tree_plot,dotplot)
#change order of taxa in moon plot
ordered_tips = full_tree$tip.label[full_tree$tip.label %in% as.character(captive_clades$cladeName)]
captive_clades_long$cladeName = factor(captive_clades_long$cladeName,levels = ordered_tips)
(dotplot <- ggplot(captive_clades_long)  +
    geom_point(aes(x=host_cat,y=cladeName,color = host_cat,size=percent_samples))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank()) +
    scale_color_manual(values=moon_color_scale) +
    scale_fill_manual(values=moon_color_scale))
dotplot %>% insert_left(tree_plot) 

captive_clades$cladeName = factor(captive_clades$cladeName,levels = ordered_tips)

taxonomy_dotplot <- ggplot(captive_clades) +
  geom_point(aes(x=0,y=cladeName,color = cladeTax))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank()) 
taxonomy_dotplot 

#add HR_type
HR_type_dotplot <- ggplot(captive_clades) +
  geom_point(aes(x=0,y=cladeName,color = HR_type))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank()) +
  scale_color_manual(values=c("dodgerblue3","darkorange2","chocolate4",
                              "chocolate4","chocolate4","purple")) 
HR_type_dotplot



all_dotplots <- plot_grid(HR_type_dotplot,taxonomy_dotplot,dotplot,nrow = 1)
ggsave('all_dotplots.pdf',height=10,width=18)

tree_plot
ggsave(tree_plot,filename='tree_plot.pdf',width=5,height=10)
