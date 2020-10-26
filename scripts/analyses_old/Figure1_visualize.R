library(tidyverse)
library(cowplot)

setwd('/Volumes/AHN/captive_ape_microbiome')

table_file <- 'results/16s/analyses/figures/16S_captive_Figure1A_data.txt'
Figure1A_data = read.table(table_file,header=TRUE,sep='\t')

(Figure1A <- ggplot(Figure1A_data, aes(numEnclosure,fill=multi_site_sp)) +
  geom_bar() + 
  ylab('# of 16S ASVs')+
  xlab('# of enclosures')+
  theme_bw()+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a'))+
  scale_x_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10)))

table_file <- 'results/16s/analyses/figures/16S_captive_Figure1B_data.txt'
Figure1B_data = read.table(table_file,header=TRUE,sep='\t')

(Figure1B <- ggplot(Figure1B_data, aes(x=sp_site_comp, y=prop_shared,fill=sp_site_comp)) + 
    geom_violin()+
    theme_bw()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared ASVs'))#

#
mean(Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='diff_spec_diff_site'])
mean(Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='diff_spec_same_site'])
mean(Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='same_spec_diff_site'])
mean(Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='same_spec_same_site'])

library(RVAideMemoire)

perm.t.test(Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='diff_spec_diff_site'],
            Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='same_spec_diff_site'])
perm.t.test(Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='diff_spec_diff_site'],
            Figure1B_data$prop_shared[Figure1B_data$sp_site_comp=='diff_spec_same_site'])

table_file <- 'results/16s/analyses/figures/16S_captive_Figure1C_data.txt'
Figure1C_data = read.table(table_file,header=TRUE,sep='\t')

Figure1C_data <-Figure1C_data %>% filter(sp_site_comp=='diff_spec_diff_site') 

(Figure1C <- ggplot(Figure1C_data, aes(x=order, y=prop_shared,fill=order)) + 
    geom_boxplot(width = 0.3,
                 position = position_dodge(width = 1),
                 outlier.shape = NA)+
    theme_bw()+
    ylab('Proportion of shared ASVs')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

alldf <- Figure1C_data %>% filter(sp_site_comp=='diff_spec_diff_site' & order=='all') 
for (x in unique(Figure1C_data$order)) {
  print(x)
  gen <- Figure1C_data %>% filter(sp_site_comp=='diff_spec_diff_site' & order==x) 
  print(perm.t.test(alldf$prop_shared,gen$prop_shared))
}

Figure1 <- plot_grid(Figure1A, Figure1B, Figure1C, nrow = 1,labels = "AUTO")
ggsave(height = 3,width=15,'results/16s/analyses/figures/Fig1.pdf')

#supplemental figures
FigureS1_sharing_bn_enclosures <- Figure1B_data %>% filter(sp_site_comp=='diff_spec_same_site') 
(ggplot(FigureSX_sharing_bn_enclosures, aes(x=full_desc_comp, y=prop_shared,fill=full_desc_comp)) + 
    geom_violin()+
    scale_fill_manual(values=c('orange2','orange2','orange2','green3',
                               'skyblue','skyblue','skyblue'
                               ))+
    theme_bw()+
    ylab('Proportion of shared ASVs')+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
ggsave('results/16s/analyses/figures/FigS1.pdf')


FigureS2_Figure1A_facet_by_HRcat <- ggplot(Figure1A_data, aes(numEnclosure,fill=multi_site_sp)) +
    geom_bar() + 
    ylab('# of 16S ASVs')+
    xlab('# of enclosures')+
    theme_bw()+
    scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a'))+
    scale_x_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10))+
    facet_wrap(~HR_cat)
FigureS2_Figure1A_facet_by_HRcat
ggsave('results/16s/analyses/figures/FigS2.pdf')

table_file <- 'results/16s/analyses/figures/16S_captive_Figure1C_data.txt'
Figure1C_data = read.table(table_file,header=TRUE,sep='\t')
Bacteroidales_data <- Figure1C_data %>% filter(order=='Bacteroidales')
(FigureS3_Bacteroidales <- ggplot(Bacteroidales_data, aes(x=sp_site_comp, y=prop_shared,fill=sp_site_comp)) + 
    geom_violin()+
    theme_bw()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared ASVs'))#
perm.t.test(Bacteroidales_data$prop_shared[Bacteroidales_data$sp_site_comp=='diff_spec_diff_site'],
            Bacteroidales_data$prop_shared[Bacteroidales_data$sp_site_comp=='same_spec_diff_site'])
perm.t.test(Bacteroidales_data$prop_shared[Bacteroidales_data$sp_site_comp=='diff_spec_diff_site'],
            Bacteroidales_data$prop_shared[Bacteroidales_data$sp_site_comp=='diff_spec_same_site'])

table_file <- 'results/16s/analyses/figures/16S_captive_FigureS3_Prevotella_data.txt'
Prevotella_data = read.table(table_file,header=TRUE,sep='\t')
(FigureS3_Prevotella <- ggplot(Prevotella_data, aes(x=sp_site_comp, y=prop_shared,fill=sp_site_comp)) + 
    geom_violin()+
    theme_bw()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared ASVs'))#
perm.t.test(Prevotella_data$prop_shared[Prevotella_data$sp_site_comp=='diff_spec_diff_site'],
            Prevotella_data$prop_shared[Prevotella_data$sp_site_comp=='same_spec_diff_site'])
perm.t.test(Prevotella_data$prop_shared[Prevotella_data$sp_site_comp=='diff_spec_diff_site'],
            Prevotella_data$prop_shared[Prevotella_data$sp_site_comp=='diff_spec_same_site'])
FigureS3 <- plot_grid(FigureS3_Bacteroidales, FigureS3_Prevotella, nrow = 1,labels = "AUTO")
ggsave(width=12,height=5,'results/16s/analyses/figures/FigS3.pdf')
