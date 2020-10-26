library("phyloseq")
library("tidyverse")
library("picante")
library("PMCMRplus")
library("cowplot")
library("reshape2")

theme_set(theme_bw())

ROOTDIR = '/Volumes/AHN/captive_ape_microbiome/'
setwd(ROOTDIR)
#INFILE
(physeq16s <- readRDS('results/16s/inputs/phyloseq_rare10000.rds'))
metadata = data.frame(sample_data(physeq16s)) %>% select(dataset,genus_sp,subspecies,common_name,
                                                         site,site_code,country,
                                                         captivity_status,captivity_status_ext,Description)
#OUTDIR
table_outdir <- 'results/16s/analyses/tables/'
dir.create(table_outdir,recursive=TRUE)
figure_outdir <-'results/16s/analyses/figures/'
dir.create(figure_outdir,recursive=TRUE)

#ALPHA DIVERSITY
calc_alpha_div <- function(physeq){
  alpha = estimate_richness(physeq, measures =  c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
  pd = picante::pd(samp = t(otu_table(physeq)), tree = phy_tree(physeq), include.root = F)
  metadata = data.frame(sample_data(physeq)) %>% select(dataset,genus_sp,subspecies,common_name,
                                                           site,site_code,country,
                                                           captivity_status,captivity_status_ext,Description)
  alphadiv <- cbind(metadata,alpha,pd) 
  alphadiv <- data.frame(alphadiv) %>% rownames_to_column(var='SampleID')
  
  #reorder captivity status for plotting
  alphadiv$captivity_status <- as.vector(recode(alphadiv$captivity_status, wild = 'wild ape', 
                                                captive = 'captive ape',
                                                industrialized_human = "industrialized human",
                                                non_industrialized_human="non-industrialized human"))
  new_order = c("wild ape","captive ape","non-industrialized human","industrialized human")
  alphadiv$captivity_status <- factor(alphadiv$captivity_status, levels = new_order)
  
  return(alphadiv)
}

#color vectors
ape_palette = c("darkorange1",'deepskyblue',"goldenrod1",'mediumorchid4')

#Visualize - alpha div
plot_alpha_figure <- function(alpha_div_table, metric) {
  alpha_plot <- ggplot(alpha_div_table, aes_string(x = "captivity_status", y = metric, fill="captivity_status"))  + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none") +
  ylab(paste0(metric," Diversity")) +
  scale_fill_manual(values = ape_palette)
  return(alpha_plot)
  }

#ALPHADIV - ALL TAXA 
alphadiv_alltaxa <- calc_alpha_div(physeq16s) 
write.table(alphadiv_alltaxa,file=file.path(table_outdir,'Table_alpha_div_16S.txt'),
            sep='\t', quote=F, row.names = F)
#PD
alphadiv_alltaxa_PD_plot <- plot_alpha_figure(alphadiv_alltaxa,'PD')
ggsave(alphadiv_alltaxa_PD_plot,file = file.path(figure_outdir,'Fig_alphadiv_16S_PD.pdf'))
kruskalTest(PD ~ captivity_status, alphadiv_alltaxa, dist="KruskalWallis")
res <- kwAllPairsDunnTest(PD ~ captivity_status, alphadiv_alltaxa,
            p.adjust.method = "fdr")
write.table(res$p.value,file.path(table_outdir,'alphadiv_16s_PD_alltaxa_dunnallpairs.txt'), 
            sep = '\t', quote = F)

#OBSERVED
alphadiv_alltaxa_Observed_plot <- plot_alpha_figure(alphadiv_alltaxa,'Observed')
ggsave(alphadiv_alltaxa_Observed_plot, file = file.path(figure_outdir,'Fig_alphadiv_16S_Observed.pdf'))
kruskalTest(PD ~ captivity_status, alphadiv_alltaxa, dist="KruskalWallis")
res <- kwAllPairsDunnTest(PD ~ captivity_status, alphadiv_alltaxa,
                          p.adjust.method = "fdr")
write.table(res$p.value,file.path(table_outdir,'alphadiv_16s_Observed_alltaxa_dunnallpairs.txt'), 
            sep = '\t', quote = F)

#select Phylum physeqs
alphdiv_by_phylum <- function(phylum_physeq16s,phylum){
  phylum_alphadiv <- calc_alpha_div(phylum_physeq16s) 
  phylum_alphadiv_Observed_plot <- plot_alpha_figure(phylum_alphadiv,'Observed')
  kwAllPairsDunnTest(Observed ~ captivity_status, phylum_alphadiv,
                     p.adjust.method = "fdr")
  res <- kwAllPairsDunnTest(PD ~ captivity_status, phylum_alphadiv,
                            p.adjust.method = "fdr")
  write.table(res$p.value,file.path(table_outdir,paste0('alphadiv_16s_Observed_',phylum,'_dunnallpairs.txt')), 
              sep = '\t', quote = F)
  return(phylum_alphadiv_Observed_plot)
  }

phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Actinobacteria')
Actinobacteria_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
(Actinobacteria_alpha_plot <- alphdiv_by_phylum(Actinobacteria_physeq16s,'Actinobacteria') + 
    ylab("# of Actinobacteria ASVs") + 
    theme(axis.text.x=element_blank()))
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Bacteroidetes')
Bacteroidetes_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
(Bacteroidetes_alpha_plot <- alphdiv_by_phylum(Bacteroidetes_physeq16s,'Bacteroidetes') + 
    ylab("# of Bacteroidetes ASVs")+ 
    theme(axis.text.x=element_blank()))
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Firmicutes') 
Firmicutes_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
(Firmicutes_alpha_plot <- alphdiv_by_phylum(Firmicutes_physeq16s,'Firmicutes') + 
    ylab("# of Firmicutes ASVs") + 
    theme(axis.text.x=element_blank()))
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Proteobacteria')
Proteobacteria_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
(Proteobacteria_alpha_plot <- alphdiv_by_phylum(Proteobacteria_physeq16s,'Proteobacteria') + 
    ylab("# of Proteobacteria ASVs") + 
    theme(axis.text.x=element_blank()))

Other_physeq16s <- subset_taxa(physeq16s,Phylum!="Actinobacteria"&
                              Phylum!="Bacteroidetes"&
                              Phylum!="Firmicutes"&
                              Phylum!="Proteobacteria")
Other_physeq16s <- prune_samples(sample_sums(Other_physeq16s)>0,Other_physeq16s)
(Other_alpha_plot <- alphdiv_by_phylum(Other_physeq16s,'Other') + 
    ylab("# of Other Phyla ASVs") + 
  theme(axis.text.x=element_blank()))

ALL_alpha_plot <- alphadiv_alltaxa_Observed_plot + 
  ylab("# of ASVs in All Phyla") + 
  theme(axis.text.x=element_blank())

alpha_by_phyla <- plot_grid(Actinobacteria_alpha_plot,Bacteroidetes_alpha_plot,
          Firmicutes_alpha_plot,Proteobacteria_alpha_plot,
          Other_alpha_plot,ALL_alpha_plot, ncol=2, labels="AUTO") 

ggsave(alpha_by_phyla,
       file = file.path(figure_outdir,'Fig_alpha_diversity_by_phyla.pdf'),
       width=6)

#BETA DIVERSITY
#Unweighted unifrac
description_color_scale <- levels(as.factor(metadata$Description))
(description_color_scale <- as.vector(recode(description_color_scale, wild_gorilla = "darkgreen", 
                                         wild_chimp = "darkorange2",
                                         wild_bonobo = "red2", 
                                         industrialized_human = "blue",
                                         non_industrialized_human="skyblue2",
                                         captive_gorilla = "darkolivegreen3", 
                                         captive_chimp = "tan1",
                                         captive_bonobo = "indianred2", 
                                         captive_orangutan="plum3")))

physeq_pcoa_unweighted_unifrac <- ordinate(
  physeq = physeq16s,
  method = "PCoA", 
  distance = "unweigthed_unifrac"
) 
plot_ordination(
  physeq = physeq16s,
  ordination = physeq_pcoa_unweighted_unifrac,
  color = "captivity_status") + 
  scale_color_manual(values = c('darkorange2',"#1F78B4","#A6CEE3",'orangered3')) +
  theme_bw()
pcoa_unweighted_unifrac <- plot_ordination(
  physeq = physeq16s,
  ordination = physeq_pcoa_unweighted_unifrac,
  color = "Description") +
  scale_color_manual(values = description_color_scale)

physeq_pcoa_weighted_unifrac <- ordinate(
  physeq = physeq16s,
  method = "PCoA", 
  distance = "weigthed_unifrac"
) 
plot_ordination(
  physeq = physeq16s,
  ordination = physeq_pcoa_weighted_unifrac,
  color = "captivity_status") + 
  scale_color_manual(values = c('darkorange2',"#1F78B4","#A6CEE3",'orangered3')) +
  theme_bw()
pcoa_weighted_unifrac <- plot_ordination(
  physeq = physeq16s,
  ordination = physeq_pcoa_weighted_unifrac,
  color = "Description") +
  scale_color_manual(values = description_color_scale)
pcoa_weighted_unifrac
pcoa_both <- plot_grid(pcoa_unweighted_unifrac,pcoa_weighted_unifrac, labels = c('A', 'B'), ncol=1)
ggsave(pcoa_both,
       file = file.path(figure_outdir,'Fig_PCoA_unweighted_weighted_unifrac.pdf'),
       width=6)
pcoa_both
#BETADISPER and PERMANOVA - Captivity status all groups
run_betadisper_permanova <- function(physeq_obj,dist_metric){
  physeq_dist <- phyloseq::distance(physeq_obj,dist_metric)
  metadata = as.data.frame(as.matrix(sample_data(physeq_obj))) %>% 
    rownames_to_column(var="X.SampleID") 
  
  print('testing homogeneity of sample groups')
  beta <- betadisper(physeq_dist, metadata$captivity_status)
  y <- permutest(beta)$tab
  rownames(y) <- paste0(dist_metric,'.',rownames(y)) 
  y = y %>% rownames_to_column(var='stat')
  colnames(y) <- paste0('betadisper_',colnames(y))
  y = rbind(colnames(y),y)
  colnames(y) <- NULL
  colnames(y) <- c('V1','V2','V3','V4','V5','V6','V7')
  print('permanova')
  z <- adonis(physeq_dist ~ captivity_status, data = metadata)
  z <- as.data.frame(z$aov.tab)
  rownames(z) <- paste0(dist_metric,'.',rownames(z)) 
  z = z %>% rownames_to_column(var='stat')
  colnames(z) <- paste0('permanova_',colnames(z))
  z = rbind(colnames(z),z)
  colnames(z) <- c('V1','V2','V3','V4','V5','V6','V7')
  return(rbind(y,z))
}

jaccard_res <- run_betadisper_permanova(physeq16s,'jaccard')
write.table(jaccard_res,'results/16s/analyses/tables/betadisper_permanova_jaccard.txt',sep='\t',quote=F,row.names = F,col.names = F)
braycurtis_res <- run_betadisper_permanova(physeq16s,'bray')
write.table(braycurtis_res,'results/16s/analyses/tables/betadisper_permanova_braycurtis.txt',sep='\t',quote=F,row.names = F,col.names = F)
unweighted_unifrac_res <- run_betadisper_permanova(physeq16s,'unweighted_unifrac')
write.table(unweighted_unifrac_res,'results/16s/analyses/tables/betadisper_permanova_unweighted_unifrac.txt',sep='\t',quote=F,row.names = F,col.names = F)
weighted_unifrac_res <- run_betadisper_permanova(physeq16s,'weighted_unifrac')
write.table(weighted_unifrac_res,'results/16s/analyses/tables/betadisper_permanova_weighted_unifrac.txt',sep='\t',quote=F,row.names = F,col.names = F)

#PAIRWISE BETADISPER and PERMANOVA 
#pairwise permanova 
cp_nih <- subset_samples(physeq16s,captivity_status =='captive'|captivity_status =='non_industrialized_human')
cp_ih <- subset_samples(physeq16s,captivity_status =='captive'|captivity_status =='industrialized_human')
cp_wd <- subset_samples(physeq16s,captivity_status =='captive'|captivity_status =='wild')
nih_ih <- subset_samples(physeq16s,captivity_status =='non_industrialized_human'|captivity_status =='industrialized_human')
nih_wd <- subset_samples(physeq16s,captivity_status =='non_industrialized_human'|captivity_status =='wild')
ih_wd <- subset_samples(physeq16s,captivity_status =='industrialized_human'|captivity_status =='wild')

jaccard_res <- lapply(c(cp_nih,cp_ih,cp_wd,nih_ih,nih_wd,ih_wd),function(x){run_betadisper_permanova(x,'jaccard')})
names(jaccard_res) <- c('cp_nih','cp_ih','cp_wd','nih_ih','nih_wd','ih_wd')
jaccard_df <- as.data.frame(do.call(rbind, jaccard_res)) %>% rownames_to_column(var='comparison')
write.table(jaccard_df,file=file.path(table_outdir,'pw_betadisper_permanova_jaccard.txt'),sep='\t',row.names = F,col.names=F)

bray_res <- lapply(c(cp_nih,cp_ih,cp_wd,nih_ih,nih_wd,ih_wd),function(x){run_betadisper_permanova(x,'bray')})
names(bray_res) <- c('cp_nih','cp_ih','cp_wd','nih_ih','nih_wd','ih_wd')
bray_df <- as.data.frame(do.call(rbind, bray_res)) %>% rownames_to_column(var='comparison')
write.table(bray_df,file=file.path(table_outdir,'pw_betadisper_permanova_braycurtis.txt'),sep='\t',row.names = F,col.names=F)

unifrac_res <- lapply(c(cp_nih,cp_ih,cp_wd,nih_ih,nih_wd,ih_wd),function(x){run_betadisper_permanova(x,'unifrac')})
names(unifrac_res) <- c('cp_nih','cp_ih','cp_wd','nih_ih','nih_wd','ih_wd')
unifrac_df <- as.data.frame(do.call(rbind, unifrac_res)) %>% rownames_to_column(var='comparison')
write.table(unifrac_df,file=file.path(table_outdir,'pw_betadisper_permanova_unweighted_unifrac.txt'),sep='\t',row.names = F,col.names=F)

wunifrac_res <- lapply(c(cp_nih,cp_ih,cp_wd,nih_ih,nih_wd,ih_wd),function(x){run_betadisper_permanova(x,'wunifrac')})
names(wunifrac_res) <- c('cp_nih','cp_ih','cp_wd','nih_ih','nih_wd','ih_wd')
wunifrac_df <- as.data.frame(do.call(rbind, wunifrac_res)) %>% rownames_to_column(var='comparison')
write.table(wunifrac_df,file=file.path(table_outdir,'pw_betadisper_permanova_weighted_unifrac.txt'),sep='\t',row.names = F,col.names=F)

#DIFFERENTIAL ABUNDANCE
collaspe_physeq <- function(physeq_obj) {
  
  physeq_obj_genus <- physeq_obj %>%  #collapse physeq to Genus level for correlations
    tax_glom("Genus") %>% 
    transform_sample_counts(function(x) x / sum(x)) #transform counts to relative abundance
  
  metadata = as.data.frame(as.matrix(sample_data(physeq_obj))) %>%  #extract metadata
    rownames_to_column(var="X.SampleID") 
  
  OTU <- data.frame(otu_table(physeq_obj_genus)) %>% #extract OTU table 
    rownames_to_column(var='ASV')
  
  TAX <- data.frame(as.matrix(tax_table(physeq_obj_genus)[, c("Phylum","Class","Order","Family","Genus")])) %>% #extract TAX table
    rownames_to_column(var='ASV') %>% 
    unite('Phylum_Genus',Phylum:Genus,sep='_')
  
  OTU_TAX <- TAX %>% 
    left_join(OTU,on = 'ASV') %>% #combine OTU and TAX table
    select(-ASV) %>%  #remove ASV column
    column_to_rownames(var = 'Phylum_Genus') %>%  #more taxonomic ranks to rownames
    t() %>% #transpose table
    as.data.frame() %>%
    rownames_to_column(var='X.SampleID') %>% #melt all genus into one column to allow for faceting
    reshape2::melt(id='X.SampleID') %>% 
    separate(variable,sep='_',into=c('Phylum','Class','Order','Family','Genus'),extra='merge')  %>% 
    dplyr::rename(relative_abundance = value) 
  print(head(OTU_TAX))
  
  OTU_TAX_META <- metadata %>% #add some metadata to sample info
    select(X.SampleID,captivity_status,site,Description) %>% 
    as.data.frame() %>% 
    left_join(OTU_TAX,by='X.SampleID')
  
  OTU_TAX_META$relative_abundance <- as.numeric(OTU_TAX_META$relative_abundance)
  
  return(OTU_TAX_META)
}

physeq_gen_abund <- collaspe_physeq(physeq16s)
physeq_gen_abund $captivity_status <- as.factor(physeq_gen_abund$captivity_status)

#calculate pairwise differences among groups using kruskal wallis
get_kw <- function(a_genus) {
  df <- physeq_gen_abund %>% filter(Genus == a_genus)
  kw.res <- kwAllPairsDunnTest(relative_abundance ~ captivity_status, df,
                               p.adjust.method = "none")
  linear.res  <- data.frame(melt(as.matrix(kw.res$p.value)))
  linear.res$value  <- str_replace(as.character(linear.res$value),'NaN','1')
  linear.res  <- linear.res %>% dplyr::filter(!is.na(value)) %>% unite(comparison, Var1:Var2, sep = 'vs.')
  linear.res$Genus <- a_genus
  return(linear.res)
}
get_kw("Bacteroides")

#calculate for all genera
df = tibble()
for (genus in unique(physeq_gen_abund$Genus)) {
  print(genus)
  res <- get_kw(genus)
  df = df %>% bind_rows(res)
}

#adjust pvalues for multiple comparisons
df$p.adj <- round(p.adjust(df$value,method = "fdr",n = length(df$value)),5)
length(unique(df$Genus))

#reshape dataframe so each row is one genus
df <- df %>% select(Genus,comparison,p.adj) %>% spread(comparison,p.adj)

#calc mean of genera by captivity status 
physeq_gen_abund_summary = physeq_gen_abund %>% 
  group_by(Genus,captivity_status) %>% 
  dplyr::summarise(mean_ra = mean(relative_abundance))
head(physeq_gen_abund_summary)

#determine the mean abundance of each genus
get_mean<- function(a_genus,group){
  select_row <- physeq_gen_abund_summary %>% filter(Genus == a_genus & captivity_status == group) 
  mean_ra <- select_row$mean_ra
  return(mean_ra)
}
get_mean('Bacteroides','wild')

#add mean abundance of each captivity status group
df$wild_mean <- mapply(get_mean, df$Genus, 'wild')
df$captive_mean <- mapply(get_mean, df$Genus, 'captive')
df$nih_mean <- mapply(get_mean, df$Genus, 'non_industrialized_human')
df$ih_mean <- mapply(get_mean, df$Genus, 'industrialized_human')

#add taxonomy 
df <- physeq_gen_abund %>% 
  select(Phylum,Class,Order,Family,Genus) %>% 
  unique() %>% 
  left_join(df,by='Genus')
dim(df) 
write.table(df,file=file.path(table_outdir,'Table_differential_abundance_genera_16S.txt'),sep='\t',row.names = F)


#genera depleted in captive apes relative to wild apes
df %>% 
  filter(wildvs.captive<.05) %>% 
  filter(wild_mean>captive_mean) %>% 
  nrow()
#genera depleted in captive apes and humans relative to wild apes
df %>% 
  filter(wildvs.captive<.05 & wildvs.non_industrialized_human<.05 & wildvs.industrialized_human<.05) %>% 
  filter(wild_mean>captive_mean & wild_mean>nih_mean & wild_mean>ih_mean) %>% 
  nrow()

#genera enriched in captive apes relative to wild apes
enriched_cp <- df %>% 
  filter(wildvs.captive<.05) %>% 
  filter(wild_mean<captive_mean) %>% 
  nrow()
#genera enriched in captive apes and humans relative to wild apes
enriched_cp_nih_ih <- df %>% 
  filter(wildvs.captive<.05 & wildvs.non_industrialized_human<.05 & wildvs.industrialized_human<.05) %>% 
  filter(wild_mean<captive_mean & wild_mean<nih_mean & wild_mean<ih_mean)  
enriched_cp_nih <- df %>% 
  filter(wildvs.captive<.05 & industrialized_humanvs.captive<.05 & non_industrialized_humanvs.industrialized_human<.05 & wildvs.industrialized_human<.05) %>% 
  filter(wild_mean<captive_mean & ih_mean<captive_mean & ih_mean<nih_mean & wild_mean<nih_mean)

#genera that reach 5% abundance in at least one group
df_over05 <-  df %>% filter(wild_mean>.05|captive_mean>.05|nih_mean>.05|ih_mean>.05)
select_gen <- physeq_gen_abund %>% filter(Genus %in% df_over05$Genus) 
select_gen$captivity_status <-ordered(select_gen$captivity_status, levels = c('wild','captive','non_industrialized_human','industrialized_human'))
select_gen$Genus <-ordered(select_gen$Genus, levels = c("Treponema_2","Prevotella_9","Ruminococcaceae_UCG-005",
                                                        "Bacteroides","Blautia","Faecalibacterium",
                                                        "Acinetobacter","Coriobacteriaceae_UCG-003","Flexilinea"))
ape_palette = c("darkorange1",'deepskyblue',"goldenrod1",'mediumorchid4')
select_gen$captivity_status <- as.vector(recode(select_gen$captivity_status, wild = 'wild ape', 
                                                       captive = 'captive ape',
                                                       industrialized_human = "industrialized human",
                                                       non_industrialized_human="non-industrialized human"))
new_order = c("wild ape","captive ape","non-industrialized human","industrialized human")
select_gen$captivity_status <- factor(select_gen$captivity_status, levels = new_order)

relabundance_plot <- ggplot(select_gen)  + 
  geom_boxplot(aes(x = captivity_status, y = relative_abundance, fill=captivity_status),outlier.shape = NA) + 
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = ape_palette) +
  facet_wrap("Genus",scale="free",nrow=3) +
  theme(axis.title.x = element_blank()) +
  ylab("Relative abundance") 
ggsave(relabundance_plot,file = file.path(figure_outdir,'Fig_relabundance_common_genera_16S.pdf'),width=9)
