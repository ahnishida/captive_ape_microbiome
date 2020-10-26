library(tidyverse)
library(ggplot2)
library(cowplot)
library(RVAideMemoire)

setwd('/Volumes/AHN/captive_ape_microbiome')
#table_file <- 'results/gyrb/analyses/codiv_moeller_ASVs/numEnclosures_table.txt'
table_file <- 'results/16s/analyses/tables/16s_prop_shared_ASVs_all_gen_table.txt'

data = read.table(table_file,header=TRUE,sep='\t')
head(data)
#p1_data <- p1_data  %>% filter(ind1_ASVs > 0)
head(data)
(p1 <- ggplot(data, aes(x=sp_site_comp, y=prop_shared_ASVs,fill=sp_site_comp)) + 
  geom_violin()+
  theme_bw()+
  scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
  ylab('Proportion of shared ASVs'))#+
 # theme(axis.title.x=element_blank(),
 #         axis.text.x=element_blank(),
 #         axis.ticks.x=element_blank()))

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
            data$prop_shared_ASVs[data$sp_site_comp=='same_spec_diff_site'])

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
            data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_same_site'])

mean(data$prop_shared_ASVs[data$sp_site_comp!='same_spec_same_site'])
mean(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'])
mean(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_same_site'])
mean(data$prop_shared_ASVs[data$sp_site_comp=='same_spec_diff_site'])
head()

subset <- data %>% filter(sp_site_comp=='diff_spec_same_site') 
subset$full_desc_comp <-  recode(subset$full_desc_comp,
  captive_orangutan_COMZ_vs_captive_gorilla_COMZ = 'COMZ_gorilla_vs_orangutan',
  captive_gorilla_COMZ_vs_captive_orangutan_COMZ= 'COMZ_gorilla_vs_orangutan',
  captive_gorilla_HOUZ_vs_captive_orangutan_HOUZ = 'HOUZ_gorilla_vs_orangutan',
  captive_chimp_HOUZ_vs_captive_orangutan_HOUZ = 'HOUZ_chimp_vs_orangutan',
  captive_chimp_HOUZ_vs_captive_gorilla_HOUZ   = 'HOUZ_chimp_vs_gorilla',
  captive_gorilla_COLZ_vs_captive_orangutan_COLZ = 'COLZ_gorilla_vs_orangutan',
  captive_bonobo_COLZ_vs_captive_orangutan_COLZ  = 'COLZ_bonobo_vs_orangutan',
  captive_bonobo_COLZ_vs_captive_gorilla_COLZ = 'COLZ_bonobo_vs_gorilla'
  )
subset$full_desc_comp <- droplevels(subset$full_desc_comp) 
subset$full_desc_comp <- factor(subset$full_desc_comp, levels = c("COLZ_bonobo_vs_gorilla",
                                                                  "COLZ_bonobo_vs_orangutan", 
                                                                  "COLZ_gorilla_vs_orangutan", 
                                                                  "HOUZ_chimp_vs_gorilla", 
                                                                  "HOUZ_chimp_vs_orangutan",
                                                                  "HOUZ_gorilla_vs_orangutan",
                                                                  "COMZ_gorilla_vs_orangutan"))
(ggplot(subset, aes(x=full_desc_comp, y=prop_shared_ASVs,fill=full_desc_comp)) + 
    geom_violin()+
    scale_fill_manual(values=c('orange2','orange2','orange2',
                               'skyblue','skyblue','skyblue',
                               'green3'))+
    theme_bw()+
    ylab('Proportion of shared ASVs')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
ggsave(height = 6,width=6,'results/16s/analyses/figures/S2.pdf')

table_file <- 'results/16s/analyses/tables/16s_prop_shared_ASVs_Bacteroidales_table.txt'
data = read.table(table_file,header=TRUE,sep='\t')
(p2 <- ggplot(data, aes(x=sp_site_comp, y=prop_shared_ASVs,fill=sp_site_comp)) + 
    geom_violin()+
    theme_bw()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared Bacteroidales ASVs')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
            data$prop_shared_ASVs[data$sp_site_comp=='same_spec_diff_site'])

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
            data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_same_site'])

table_file <- 'results/16s/analyses/tables/16s_prop_shared_ASVs_Prevotella_table.txt'
data = read.table(table_file,header=TRUE,sep='\t')
(p3 <- ggplot(data, aes(x=sp_site_comp, y=prop_shared_ASVs,fill=sp_site_comp)) + 
    geom_violin()+
    theme_bw()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared Prevotella ASVs')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
            data$prop_shared_ASVs[data$sp_site_comp=='same_spec_diff_site'])

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
            data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_same_site'])

table_file <- 'results/gyrb/analyses/codiv_moeller_ASVs/gyrb_prop_shared_ASVs_Prevotella_table.txt'
data = read.table(table_file,header=TRUE,sep='\t')
(p4 <- ggplot(data, aes(x=sp_site_comp, y=prop_shared_ASVs,fill=sp_site_comp)) + 
    geom_violin()+
    theme_light()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared ASVs')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

table_file <- 'results/16s/analyses/tables/16s_prop_shared_ASVs_top10_table.txt'
data = read.table(table_file,header=TRUE,sep='\t')
subset <- data %>% filter(sp_site_comp=='diff_spec_diff_site') 
head(subset)
dodge <- position_dodge(width = 1)

(p5 <- ggplot(subset, aes(x=order, y=prop_shared_ASVs,fill=order)) + 
    geom_boxplot(width = 0.3,
                 position = dodge,
                 outlier.shape = NA)+
    theme_bw()+
    #scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared ASVs')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))
unique(data$order)
alldf <- data %>% filter(sp_site_comp=='diff_spec_diff_site' & order=='all') 
gen <- data %>% filter(sp_site_comp=='diff_spec_diff_site' & order=='Bacteroidales') 
perm.t.test(alldf$prop_shared_ASVs,gen$prop_shared_ASVs)
gen <- data %>% filter(sp_site_comp=='diff_spec_diff_site' & order=='Clostridiales') 
perm.t.test(alldf$prop_shared_ASVs,gen$prop_shared_ASVs)
gen <- data %>% filter(sp_site_comp=='diff_spec_diff_site' & order=='Clostridiales') 
perm.t.test(alldf$prop_shared_ASVs,gen$prop_shared_ASVs)
for (x in unique(data$order)) {
  print(x)
  gen <- data %>% filter(sp_site_comp=='diff_spec_diff_site' & order==x) 
  print(perm.t.test(alldf$prop_shared_ASVs,gen$prop_shared_ASVs))
  }



table_file <- 'results/16s/analyses/tables/numEnclosures_table.txt'
data = read.table(table_file,header=TRUE)
(e1 <- ggplot(data, aes(fill=multi_site_sp, y=count, x=numEnclosure)) + 
  geom_bar(position="stack", stat="identity") +
    ylab('# of 16S ASVs')+
    xlab('# of enclosures')+
    theme_bw()+
    scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a'))+
    scale_x_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10))) 

data %>% select(count) %>% sum()
data %>% filter(numEnclosure >= 4) %>% select(count) %>% sum()


table_file <- 'results/16s/analyses/tables/numEnclosures_Bacteroidales_table.txt'
data = read.table(table_file,header=TRUE)
(e2 <- ggplot(data, aes(fill=multi_site_sp, y=count, x=numEnclosure)) + 
    geom_bar(position="stack", stat="identity") +
    geom_bar(position="stack", stat="identity") +
    ylab('# of Bacteroidales 16S ASVs')+
    xlab('# of enclosures')+
    theme_bw() +
    scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a'))+
    scale_x_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10))) 
data %>% select(count) %>% sum()
data %>% filter(numEnclosure >= 4) %>% select(count) %>% sum()


table_file <- 'results/16s/analyses/tables/numEnclosures_Prevotella_table.txt'
data = read.table(table_file,header=TRUE)
(e3 <- ggplot(data, aes(fill=multi_site_sp, y=count, x=numEnclosure)) + 
    geom_bar(position="stack", stat="identity") +
    geom_bar(position="stack", stat="identity") +
    ylab('# of Prevotella 16S ASVs')+
    xlab('# of enclosures')+
    theme_bw() +
    scale_fill_manual(values=c('#1b9e77','#7570b3','#e7298a'))+
    scale_x_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10))) 

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'][1:438],
            data$prop_shared_ASVs[data$sp_site_comp=='same_spec_diff_site'])

perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
          data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_same_site'])
perm.t.test(data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'],
       data$prop_shared_ASVs[data$sp_site_comp=='diff_spec_diff_site'])

plot_grid(p2, p3, e2, e3, labels = "AUTO")
ggsave(height = 6,width=10,'results/16s/analyses/figures/S3.pdf')
plot_grid(p1, e1, p5, nrow = 1,labels = "AUTO")
ggsave(height = 3,width=15,'results/16s/analyses/figures/Fig1.pdf')
