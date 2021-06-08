library("phyloseq")
library("tidyverse")
library("picante")
library("PMCMRplus")
library("cowplot")
library("ggplot2")
library("reshape2")
library("RVAideMemoire")
library("RColorBrewer")
library('rstatix')
library('broom')
set.seed(3113)
#INPUTS
physeq16s <- readRDS('results/16s/inputs/phyloseq_rare10000.rds') #read in physeq
HR_table <- read.table('results/16s/analyses/tables/16S_ASVs_summary.txt',sep='\t',header=T)

#OUTDIR
table_outdir <- 'results/16s/analyses/tables' #specify output folder for tables
dir.create(table_outdir,recursive=TRUE)
figure_outdir <-'results/16s/analyses/figures' #specify output folder for figures
dir.create(figure_outdir,recursive=TRUE)

#ADD HR INFO to PHYSEQ TAX_TABLE
new_tax_table <- data.frame(tax_table(physeq16s)) #add HRcat and HRtype to taxonomy
new_tax_table$HR_cat <- HR_table$HR_cat
new_tax_table$HR_type <- HR_table$HR_type
tax_table(physeq16s) <- as.matrix(new_tax_table)

#Make new column in metadata that is a combo of description/zoo site
metadata <- sample_data(physeq16s) #format metadata with new column
metadata$Description_country_zoo <- paste(metadata$Description,metadata$country_zoo,sep='_')
sample_data(physeq16s) <- metadata

#set colors for sample groups
description_color_scale <- levels(as.factor(metadata$Description))
description_color_scale <- as.vector(recode(description_color_scale, wild_gorilla = "darkgreen",
                                             wild_chimp = "darkorange2",
                                             wild_bonobo = "red2",
                                             industrialized_human = "blue",
                                             non_industrialized_human="skyblue2",
                                             captive_gorilla = "darkolivegreen3",
                                             captive_chimp = "tan1",
                                             captive_bonobo = "indianred2",
                                             captive_orangutan="plum3"))

#ORDINATION
print('calculating beta diversity distance metrics')
physeq16s_bray <- phyloseq::distance(physeq16s,'bray')
physeq16s_jaccard <- phyloseq::distance(physeq16s,'jaccard')
physeq16s_wunifrac <- phyloseq::distance(physeq16s,'weighted_unifrac')
physeq16s_uunifrac <- phyloseq::distance(physeq16s,'unweighted_unifrac')

ordinate_NMDS <- function(physeq,physeq_dist){
  #generate NMDS ordination plot
  ord_nMDS<- metaMDS(comm= physeq_dist,
                          k = 2,
                          maxit = 1000, #have to increase tries for unifrac distances, otherwise they don't converge
                          trymax = 500,
                          wascores = FALSE,
                          autotransform = FALSE,
                          trace = 2,
                          noshare = FALSE,
                          parallel = 4)
  ndim_stress = paste0('k = ',ord_nMDS$ndim,', stress = ',round(ord_nMDS$stress,2)) #extract stress and dimensions
  
  plot_ord <- plot_ordination(
    physeq = physeq,
    ordination =ord_nMDS,
    color = "Description") +
    scale_color_manual(values = description_color_scale)+
    theme_bw() +
    theme(legend.position="bottom") +
    annotate(geom="text", x=0, y=0,
             label=ndim_stress,
             color="red")
  return(plot_ord)
}

print('NMDS ordinations')
capture.output(plot_nmds_bray <- ordinate_NMDS(physeq16s,physeq16s_bray),file='nul') 
ggsave(plot_nmds_bray, file = file.path(figure_outdir,'Fig5_nMDS_bray.pdf'),width=6,height=5)
capture.output(plot_nmds_jaccard <- ordinate_NMDS(physeq16s,physeq16s_jaccard),file='nul')
capture.output(plot_nmds_wunifrac <- ordinate_NMDS(physeq16s,physeq16s_wunifrac),file='nul')
capture.output(plot_nmds_uunifrac <- ordinate_NMDS(physeq16s,physeq16s_uunifrac),file='nul')
plot_nmds_othermetrics <- plot_grid(plot_nmds_jaccard,
                                    plot_nmds_wunifrac,
                                    plot_nmds_uunifrac,
                                    ncol=1,
                                    labels = "AUTO")
ggsave(plot_nmds_othermetrics, file = file.path(figure_outdir,'FigS8_nMDS_othermetrics.pdf'),width=6,height=14)

#BETADISPER and PERMANOVA - Description all groups
run_betadisper_permanova_dist <- function(physeq,physeq_dist,dist_metric){
  #output vector containing betadisper and permanova based on sample groups
  metadata = as.data.frame(as.matrix(sample_data(physeq))) %>% #extract metadata
    rownames_to_column(var="X.SampleID")
  print('testing homogeneity of sample groups')
  beta <- betadisper(physeq_dist, metadata$Description) #run betadisper
  beta_tab <- permutest(beta)$tab #extract table
  print(beta_tab)
  beta_df <- beta_tab$Df[[1]]
  beta_F <- beta_tab$F[[1]]
  beta_p <- beta_tab$P[[1]]
  print('permanova')
  perm <- adonis(physeq_dist ~ Description, data = metadata) #run permanova
  perm_tab <- as.data.frame(perm$aov.tab) #extract table
  print(perm_tab)
  perm_pseudof <- perm_tab$F.Model[[1]]
  perm_R2 <- perm_tab$R2[[1]]
  perm_pvalue <- perm_tab$Pr[[1]]
  beta_perm_res <- c(dist_metric,beta_df,beta_F,beta_p,perm_pseudof,perm_R2,perm_pvalue) #combine betadisper and permanova
  names(beta_perm_res) <- c('dist_metric','Df','betadisper_F','betadisper_Pr(>F)',
                            'permanova_F.Model','permanova_R2','permanova_Pr(>F)')
  return(beta_perm_res)
}

print('BETADISPER and PERMANOVA')
print('bray-curtis')
(beta_perm_bray <- run_betadisper_permanova_dist(physeq16s,physeq16s_bray,'bray-curtis'))
print('jaccard')
(beta_perm_jaccard <- run_betadisper_permanova_dist(physeq16s,physeq16s_jaccard,'jaccard'))
print('weighted unifrac')
(beta_perm_wunifrac <- run_betadisper_permanova_dist(physeq16s,physeq16s_wunifrac,'weighted unifrac'))
print('unweighted unifrac')
(beta_perm_uunifrac <- run_betadisper_permanova_dist(physeq16s,physeq16s_uunifrac,'unweighted unifrac'))
#combine permanova from distance metrics to generate table
permanova_beta_tab <- data.frame(rbind(beta_perm_bray,beta_perm_jaccard,beta_perm_wunifrac,beta_perm_uunifrac))
permanova_beta_tab <- permanova_beta_tab %>% mutate(betadisper_F = round(as.numeric(betadisper_F),1),
                                                    permanova_F.Model = round(as.numeric(permanova_F.Model),1),
                                                    permanova_R2= round(as.numeric(permanova_R2),2))
write.table(permanova_beta_tab,file=file.path(table_outdir,'Table_betadisper_permanova_Figure5_FigureS8.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#PAIRWISE BETADISPER and PERMANOVA
#extract physeq for each of the pairwise group comparisons
Descriptions <- c("captive_orangutan","captive_gorilla","captive_chimp","captive_bonobo",
                  "wild_gorilla","wild_chimp","wild_bonobo","non_industrialized_human","industrialized_human")
pw_comps <- combn(Descriptions,2) #generate list of all combinations of 2 description
pw_comps <- as.data.frame(t(pw_comps)) #format into dataframe
colnames(pw_comps) <- c('group1','group2')
pw_df = data.frame()
print('calculate pairwise permanova distances between groups')
for (i in 1:nrow(pw_comps)){
  #I'm bad at apply functions in R, loop through combinations of sample groups
  pw <- pw_comps[i,]
  print(pw) #pairwise comparison
  physeq_pw <- subset_samples(physeq16s,Description==pw$group1|Description==pw$group2) #subset physeq to two groups
  physeq_pw_dist <- phyloseq::distance(physeq_pw,'bray') #calculate distance
  capture.output(beta_perm_res <- run_betadisper_permanova_dist(physeq_pw,physeq_pw_dist,'bray-curtis'),file='nul') #run betadisper and permanova
  res <- c(pw$group1,pw$group2,beta_perm_res) #generate new line with group names and results
  pw_df <- rbind(pw_df,res) #bind new row to dataframe
}
system('rm nul')
colnames(pw_df) <- c('group1','group2','dist_metric','Df','betadisper_F','betadisper_Pr(>F)',
                     'permanova_F.Model','permanova_R2','permanova_Pr(>F)')
pw_df <- pw_df %>% mutate(betadisper_F = round(as.numeric(betadisper_F),1),
                          permanova_F.Model = round(as.numeric(permanova_F.Model),1),
                          permanova_R2= round(as.numeric(permanova_R2),2))
write.table(pw_df,file=file.path(table_outdir,'TableS4_pairwise_betadisper_permanova_Figure5.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#PERMANOVA captive apes host filtering

#filter to just captive ape samples 
physeq16s_captive = subset_samples(physeq16s,captivity_status=='captive')
captive_metadata = as.data.frame(as.matrix(sample_data(physeq16s_captive)))
captive_metadata = captive_metadata %>% mutate(enclosure = paste(Description,site_code,sep='_'))
physeq16s_captive_bray <- phyloseq::distance(physeq16s_captive,'bray')

#betadisper
beta <- betadisper(physeq16s_captive_bray, captive_metadata$site_code) #run betadisper
(beta_tab <- permutest(beta)$tab) #
beta <- betadisper(physeq16s_captive_bray, captive_metadata$enclosure) 
(beta_tab <- permutest(beta)$tab) 
beta <- betadisper(physeq16s_captive_bray, captive_metadata$Description) 
(beta_tab <- permutest(beta)$tab) 
beta <- betadisper(physeq16s_captive_bray, captive_metadata$dataset) 
(beta_tab <- permutest(beta)$tab) 
perm <- adonis(physeq16s_captive_bray ~ site_code + enclosure + Description + dataset, data = captive_metadata)
perm

#COMPOSITION BY PHYLA
print('merging samples by description to visualize average relative abundances of phyla')
physeq16s_Description_Mean <-  merge_samples(physeq16s, "Description") #collaspe samples into groups based on species and captivity status
otu_table(physeq16s_Description_Mean) <- otu_table(physeq16s_Description_Mean)/sample_sums(physeq16s_Description_Mean) #format count to relative abundance
physeq16s_Description_Mean <- psmelt(physeq16s_Description_Mean) #collaspe to dataframe with metadata
sample_ordered <- c("wild_bonobo","wild_chimp","wild_gorilla","captive_bonobo","captive_chimp",
                    "captive_gorilla","captive_orangutan",
                    "non_industrialized_human","industrialized_human")
physeq16s_Description_Mean$Sample <- factor(physeq16s_Description_Mean$Sample, levels =sample_ordered)

print('generating list of phyla that reach mean relative abundance of 1% in one sample group')
options(dplyr.summarise.inform=F)
(Phylum_over01 <- physeq16s_Description_Mean %>% 
  group_by(Phylum) %>%
  summarise(max_abundance=max(Abundance)) %>%
  as.data.frame() %>%
  filter(max_abundance>.01) %>% 
  select(Phylum) %>% unique() %>% unlist())

print('summary of the mean relative abundance of bacterial phylum across sample groups')
(physeq16s_Description_Mean %>% 
    filter(Phylum %in% Phylum_over01) %>%
    group_by(Phylum,Sample) %>% 
    summarise(mean_abundance=sum(Abundance)) %>%
    as.data.frame() %>% 
    spread(key = Phylum, value=mean_abundance))
print('summary of the mean relative abundance of HR ASVs across sample groups')
(physeq16s_Description_Mean %>% 
  group_by(HR_type,Sample) %>% 
  summarise(mean_abundance=sum(Abundance)) %>%
  as.data.frame() %>% 
  spread(key = HR_type, value=mean_abundance))

#Frequency of conspecific ASVs 
physeq16S_melt = psmelt(physeq16s)
head(physeq16S_melt)
conspec = physeq16S_melt %>% 
  select(OTU,Description,Abundance) %>% 
  filter(Abundance>0) %>%
  group_by(OTU,Description) %>% 
  tally() %>% 
  as.data.frame() %>%
  pivot_wider(names_from = Description, values_from = n,values_fill=0) %>%
  mutate_at(vars(-OTU),as.numeric) %>%
  mutate(total = rowSums(select(conspec,-OTU)))
chimp_conspec = conspec %>%
  filter(captive_chimp+wild_chimp == total, captive_chimp>0, wild_chimp>0)
gorilla_conspec = conspec %>%
  filter(captive_gorilla+wild_gorilla== total, captive_gorilla>0, wild_gorilla>0)
bonobo_conspec = conspec %>%
  filter(captive_bonobo+wild_bonobo == total, captive_bonobo>0, wild_bonobo>0)

physeq16S_melt %>% filter(OTU %in% c(gorilla_conspec$OTU,chimp_conspec$OTU)) %>% filter(captivity_status=='captive',Abundance>0) 

#Figure 3A
(plot_phyla_abund <- physeq16s_Description_Mean %>%
  group_by(Sample,Phylum) %>% #group all ASVs within a HR_type
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  filter(Phylum %in% Phylum_over01) %>%
  ggplot(aes(fill=Phylum, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Set3") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0)))

#COMPOSITION BY PHYLA - STATS
physeq16s_stats <- physeq16s #make copy of physeq
otu_table(physeq16s_stats) <- otu_table(physeq16s_stats)/sample_sums(physeq16s_stats) #format count to relative abundance
physeq16s_melt <- psmelt(physeq16s_stats) #combines taxa table, metadata, and otu table into df

kruskal_captivitystatus <- function(TaxName,group1,group2,df) {
  #runs kruskal-wallis and retrieves group summary stats 
  group1_mean = df %>% filter(captivity_status==group1|Description==group1) %>% summarise(group1_mean = mean(Abundance)) 
  group2_mean = df %>% filter(captivity_status==group2|Description==group2) %>% summarise(group2_mean = mean(Abundance))
  group1_n = df %>% filter(captivity_status==group1|Description==group1) %>% nrow()
  group2_n = df %>% filter(captivity_status==group2|Description==group2) %>% nrow()
  if (group1_mean>0|group2_mean>0){ #check to make sure values in wild and captive apes
    kruskal <- df %>% 
        kruskal_test(Abundance ~ captivity_status) %>% 
        select(df,statistic,p) %>% 
        unlist() } else {
    kruskal <- c(NA,NA,NA)}
  df.kruskal <- unlist(c(TaxName,group1,group2,group1_n,group2_n,group1_mean,group2_mean,kruskal))
  names(df.kruskal) <- c('TaxName','group1','group2','group1_n','group2_n',
                            'group1_mean','group2_mean','df','statistic','p')
  return(df.kruskal)
}

kw_abund_Taxa <- function(TaxLevel,TaxName) {
  #runs kruskal-wallis on abundance of any bacterial taxa at any specified taxonomic level
  #runs 4 comparisons of abundance: captive vs wild apes, captive vs wild bonobo, 
  #captive vs wild chimp, captive vs wild gorilla
  
  physeq16s_Taxa <- physeq16s_melt %>%
    filter(!!sym(TaxLevel)==TaxName) %>% #subset to the just taxa of interest
    group_by(Sample,captivity_status,Description) %>%  #group all ASVs within an individual, keep group metadata
    summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
    as.data.frame()
   
  #captive vs wild all species
  cp_wd <- physeq16s_Taxa %>% filter(captivity_status=='captive'|captivity_status=='wild')
  cp_wd.kruskal <- kruskal_captivitystatus(TaxName,'captive','wild',cp_wd)
  
  #captive vs wild bonobo
  bonobo <- physeq16s_Taxa %>% filter(Description=='captive_bonobo'|Description=='wild_bonobo')
  bonobo.kruskal <- kruskal_captivitystatus(TaxName,'captive_bonobo','wild_bonobo',bonobo)
  
  #captive vs wild chimp
  chimp <- physeq16s_Taxa %>% filter(Description=='captive_chimp'|Description=='wild_chimp')
  chimp.kruskal <- kruskal_captivitystatus(TaxName,'captive_chimp','wild_chimp',chimp)
  
  #captive vs wild gorilla
  gorilla <- physeq16s_Taxa %>% filter(Description=='captive_gorilla'|Description=='wild_gorilla')
  gorilla.kruskal <- kruskal_captivitystatus(TaxName,'captive_gorilla','wild_gorilla',gorilla)
  
  res <- data.frame(rbind(cp_wd.kruskal,bonobo.kruskal,chimp.kruskal,gorilla.kruskal)) %>% 
    rownames_to_column('comparison')
  return(res)
}

print('kruskal-wallis differential abundance test')
kw_abund_Taxa('Phylum','Actinobacteria') #test on a single phylum
kw_abund_Taxa('Genus','Bifidobacterium') #test on a single genus

#run Kruskal-Wallis on all phyla over 1%
kw_abund_Phyla <- function(TaxName){
  kw_abund_Taxa('Phylum',TaxName) #helper function bc again I'm bad at apply functions
}
kw_phylum <- lapply(Phylum_over01,kw_abund_Phyla) 
kw_phylum_df <- as.data.frame(do.call(rbind,kw_phylum)) #concat into one dataframe 
kw_phylum_df$p_adj_bonferroni <- round(p.adjust(kw_phylum_df$p,method = "bonferroni"),3)
kw_phylum_df <- kw_phylum_df %>% rename(Phylum = TaxName) %>%
                 mutate(p=round(as.numeric(p),3),
                        statistic=round(as.numeric(statistic),1),
                        group1_mean=round(as.numeric(group1_mean),4),
                        group2_mean=round(as.numeric(group2_mean),4),
                        group1 = recode(group1,'captive'='captive_ape','wild'='wild_ape'),
                        group2 = recode(group2,'captive'='captive_ape','wild'='wild_ape'))
write.table(as.data.frame(kw_phylum_df),
            file=file.path(table_outdir,'Table_kruskal_wallis_phyla_Figure3A.txt'),
            sep='\t',row.names = F,col.names=T,quote=F)

#COMPOSITION BY HR_TYPE
HR_type_ordered <- c("Unique_CP","MX_human_wild_apes","MX_wild_apes",
                     "HR_human","HR_wild_gorilla","HR_wild_chimp","HR_wild_bonobo")
physeq16s_Description_Mean$HR_type <- factor(physeq16s_Description_Mean$HR_type, levels = HR_type_ordered )
HRpalette <- as.vector(recode(HR_type_ordered, HR_wild_gorilla = "darkgreen",
                               HR_wild_chimp = "darkorange2",
                               HR_wild_bonobo = "red2",
                               HR_human = "dodgerblue",
                               MX_human_wild_apes = "goldenrod",
                               MX_wild_apes = "maroon",
                               Unique_CP='purple'))
#Figure 3B
(plot_phyla_abund_HR <- physeq16s_Description_Mean %>%
  group_by(Sample,HR_type) %>% #group all ASVs within a HR_type
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=HR_type, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = HRpalette) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0)))


#COMPOSITION BY HR_TYPE - STATS
HR_type_stats <-function(physeq_melt){
  #compares relative abundance of Host-restricted and mixed ASVs across sample groups 
  df <- physeq_melt %>%
    group_by(Sample,captivity_status,Description,HR_type) %>%  #group all ASVs within an individual, keep group metadata
    summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
    as.data.frame() %>%
    mutate(captivity_status = recode(captivity_status,
                                     'non_industrialized_human'='human', 
                                     'industrialized_human'='human')) 
  #subset to human host-restricted ASVs
  HR_human_df <- df %>% filter(HR_type == 'HR_human' & captivity_status != 'wild')
  (HR_human.kruskal <- kruskal_captivitystatus('HR_human','captive','human',HR_human_df))
  #subset to wild gorilla host-restricted ASVs
  HR_wild_gorilla_df  <- df %>% filter(HR_type == 'HR_wild_gorilla') %>%
    filter(captivity_status == "captive" | Description == "wild_gorilla")
  (HR_wild_gorilla.kruskal <-  kruskal_captivitystatus('HR_wild_gorilla','captive','wild_gorilla',HR_wild_gorilla_df))
  #subset to wild chimp host-restricted ASVs
  HR_wild_chimp_df  <- df %>% filter(HR_type == 'HR_wild_chimp') %>%
    filter(captivity_status == "captive" | Description == "wild_chimp")
  (HR_wild_chimp.kruskal  <-  kruskal_captivitystatus('HR_wild_chimp','captive','wild_chimp',HR_wild_chimp_df))
  #subset to wild bonobo host-restricted ASVs
  HR_wild_bonobo_df  <- df %>% filter(HR_type == 'HR_wild_bonobo') %>%
    filter(captivity_status == "captive" | Description == "wild_bonobo")
  (HR_wild_bonobo.kruskal  <-  kruskal_captivitystatus('HR_wild_bonobo','captive','wild_bonobo',HR_wild_bonobo_df))
  #subset to mixed-host ASVs found in wild apes
  MX_wild_apes_df  <- df %>% filter(HR_type == 'MX_wild_apes' & captivity_status != 'human') 
  (MX_wild_apes.kruskal  <-  kruskal_captivitystatus('MX_wild_apes','captive','wild',MX_wild_apes_df))
  #subset to mixed-host ASVs found in wild apes and humans
  MX_human_wild_apes_df <- df %>% filter(HR_type == 'MX_human_wild_apes' & captivity_status != 'human')
  (MX_human_wild_apes.kruskal1  <-  kruskal_captivitystatus('MX_human_wild_apes','captive','wild',MX_human_wild_apes_df))
  MX_human_wild_apes_df <- df %>% filter(HR_type == 'MX_human_wild_apes' & captivity_status != 'wild')
  (MX_human_wild_apes.kruskal2  <-  kruskal_captivitystatus('MX_human_wild_apes','captive','human',MX_human_wild_apes_df))
  MX_human_wild_apes_df <- df %>% filter(HR_type == 'MX_human_wild_apes' & captivity_status != 'captive')
  (MX_human_wild_apes.kruskal3  <-  kruskal_captivitystatus('MX_human_wild_apes','human','wild',MX_human_wild_apes_df))
  #combine outputs into dataframe
  res <- data.frame(rbind(HR_human.kruskal,
                          HR_wild_gorilla.kruskal,
                          HR_wild_chimp.kruskal,
                          HR_wild_bonobo.kruskal,
                          MX_wild_apes.kruskal,
                          MX_human_wild_apes.kruskal1,
                          MX_human_wild_apes.kruskal2,
                          MX_human_wild_apes.kruskal3))
  res$p_adj_bonferroni <- round(p.adjust(res$p,method = "bonferroni"),3) #correct p-values for multiple testing
  res <- res %>% rename(HR_type = TaxName) %>%
                 mutate(p=round(as.numeric(p),3),
                          statistic=round(as.numeric(statistic),1),
                          group1_mean=round(as.numeric(group1_mean),4),
                          group2_mean=round(as.numeric(group2_mean),4),
                          group1 = recode(group1,'captive'='captive_ape','wild'='wild_ape'),
                          group2 = recode(group2,'captive'='captive_ape','wild'='wild_ape'))
  return(res)
  
}
kruskal_wallis_HR_type <- HR_type_stats(physeq16s_melt)
write.table(kruskal_wallis_HR_type,file=file.path(table_outdir,'Table_kruskal_wallis_HRtype_Figure3B.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#BACTERIAL GENERA DIFF ABUNDANCE IN WILD AND CAPTIVE APES
Genera <- unique(physeq16s_melt$Genus) #list of unique genera
Genera <- Genera[!is.na(Genera)] #remove NA from list
kw_abund_Genera <- function(TaxName){
  kw_abund_Taxa('Genus',TaxName)
}
kw_Genera_df <- lapply(Genera,kw_abund_Genera)
kw_Genera_df <- as.data.frame(do.call(rbind, kw_Genera_df)) #combine into dataframe
kw_Genera_df <- kw_Genera_df %>% rename(Genus = TaxName) %>%
                                 mutate(group1_mean = as.numeric(group1_mean), #do.call changes variable type have to apply as.numeric
                                        group2_mean = as.numeric(group2_mean),
                                        p = as.numeric(p),
                                        statistic=round(as.numeric(statistic),1),
                                        group1_mean=round(as.numeric(group1_mean),4),
                                        group2_mean=round(as.numeric(group2_mean),4),
                                        group1 = recode(group1,'captive'='captive_ape','wild'='wild_ape'),
                                        group2 = recode(group2,'captive'='captive_ape','wild'='wild_ape'))
kw_Genera_df$p_adj_bonferroni <- round(p.adjust(kw_Genera_df$p,method = "bonferroni"),3)
kw_Genera_df$p <- round(kw_Genera_df$p,3)
write.table(as.data.frame(kw_Genera_df),file=file.path(table_outdir,'TableS5_kruskal_wallis_genera.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#Identify bacterial genera that show parallel enrichment in all captive ape species
Genera_enriched <- kw_Genera_df %>% 
                              filter(group1_mean > group2_mean) %>% #captive mean > wild mean
                              filter(p_adj_bonferroni < .05) %>% #significant
                              filter(group1_mean>.01) %>%  #mean abundance must be greater than 1%
                              group_by(Genus) %>% tally() %>% filter(n==4) %>% #significant across all comparison
                              select(Genus) %>% unlist()
kw_Genera_df_enriched <- kw_Genera_df %>% filter(Genus %in% Genera_enriched)
write.table(as.data.frame(kw_Genera_df_enriched),file=file.path(table_outdir,'Table_kruskal_wallis_genera_enriched_Figure3C.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#Identify bacterial genera depleted in captive vs. wild in one comparison
Genera_depleted <- kw_Genera_df %>% 
  filter(group1_mean < group2_mean) %>% #captive mean < wild mean
  filter(p_adj_bonferroni < .01) %>%  #significant
  filter(group2_mean > .02) %>% #mean abundance must be greater than 1%
  filter(comparison != 'cp_wd.kruskal') %>% #exclude comps significant b/n all captive apes and all wild apes 
  select(Genus) %>% unique() %>% unlist()
kw_Genera_df_depleted <- kw_Genera_df %>% filter(Genus %in% Genera_depleted)
write.table(as.data.frame(kw_Genera_df_depleted),file=file.path(table_outdir,'Table_kruskal_wallis_genera_depleted_FigureS9A.txt'),sep='\t',row.names = F,col.names=T,quote=F)

print('run HR stats on genera enriched/depleted')
physeq16s_melt_Genera_enriched <- physeq16s_melt %>%
  filter(Genus %in% Genera_enriched)
(kruskal_wallis_HRtype_Genera_enriched <- HR_type_stats(physeq16s_melt_Genera_enriched))
write.table(kruskal_wallis_HRtype_Genera_enriched,file=file.path(table_outdir,'Table_kruskal_wallis_HRtype_genera_enriched_Figure3D.txt'),sep='\t',row.names = F,col.names=T,quote=F)

print('HR stats for genera depleted in captive apes')
physeq16s_melt_Genera_depleted <- physeq16s_melt %>%
  filter(Genus %in% Genera_depleted)
HR_type_stats(physeq16s_melt_Genera_depleted)
(kruskal_wallis_HRtype_Genera_depleted <- HR_type_stats(physeq16s_melt_Genera_depleted))
write.table(kruskal_wallis_HRtype_Genera_depleted,file=file.path(table_outdir,'Table_kruskal_wallis_HRtype_genera_depleted_FigureS9B.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#Figure 3C,Visualize bacterial genera enriched in captivity 
(plot_genus_enriched_abund <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_enriched) %>% #select only genera enriched in captive apes
  group_by(Sample,Genus) %>% #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=Genus, y=Abundance, x=Sample)) + #color by ASV Genus
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Dark2") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

#Figure3D: Visualize by HRtype
(plot_genus_enriched_abund_HRtype <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_enriched) %>% #select only genera enriched in captive apes
  group_by(Sample,HR_type) %>%  #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=HR_type, y=Abundance, x=Sample)) + #color by ASV HR_type
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = HRpalette) + 
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance") +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

#Figure S9A: Visualize bacterial genera depleted in captivity
colors = c(brewer.pal(name="Paired", n = 12),brewer.pal(name="Dark2", n = 8))
(plot_genus_deplete_abund <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_depleted) %>% #select Genera depleted in captive apes
  group_by(Sample,Genus) %>% #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=Genus, y=Abundance, x=Sample)) + #color by ASV genus
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

#Figure S9B: Visualize by HRtype
(plot_genus_deplete_abund_HRtype <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_depleted) %>% #select Genera depleted in captive apes
  group_by(Sample,HR_type) %>% #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=HR_type, y=Abundance, x=Sample)) +  #color by ASV HR_type
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = HRpalette) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

#Stitch together plots for Figure 3
top_row <- plot_grid(plot_phyla_abund + theme(legend.position = "none"),get_legend(plot_phyla_abund),
            plot_phyla_abund_HR + theme(legend.position = "none"),get_legend(plot_phyla_abund_HR),
            nrow=1,
            rel_widths = c(1,.75,1,.75),
            labels = c('A','','B'))
bottom_row <-plot_grid(plot_genus_enriched_abund + theme(legend.position = "none"),get_legend(plot_genus_enriched_abund),
                plot_genus_enriched_abund_HRtype + theme(legend.position = "none"),get_legend(plot_genus_enriched_abund_HRtype),
                nrow = 1,
                rel_widths = c(1,.75,1,.75),
                labels = c('C','','D'))
composition_taxa_HR <- plot_grid(top_row,bottom_row,nrow=2)
ggsave(composition_taxa_HR,file=file.path(figure_outdir,'Fig3_composition_taxa_HR.pdf'),height=11,width=12)

#Stitch together plots for Figure S9
plot_genus_deplete_abund_combined <- plot_grid(plot_genus_deplete_abund + theme(legend.position = "none"),get_legend(plot_genus_deplete_abund),
                                               plot_genus_deplete_abund_HRtype + theme(legend.position = "none"),get_legend(plot_genus_deplete_abund_HRtype),
                                               rel_widths = c(1,.5,1,.5))
ggsave(plot_genus_deplete_abund_combined,file=file.path(figure_outdir,'FigS9_genera_depleted_captive.pdf'),height=13,width=8)


#Figure S10: COMPOSITION BARPLOT BY SITE
physeq16s_DescriptionSite_Mean <-  merge_samples(physeq16s, "Description_country_zoo") #collaspe samples into groups based on species and captivity status
otu_table(physeq16s_DescriptionSite_Mean) <- otu_table(physeq16s_DescriptionSite_Mean)/sample_sums(physeq16s_DescriptionSite_Mean) #format count to relative abundance
physeq16s_DescriptionSite_Mean <- psmelt(physeq16s_DescriptionSite_Mean)
sample_ordered <- c("wild_bonobo","wild_chimp","wild_gorilla","captive_bonobo","captive_chimp",
                    "captive_gorilla","captive_orangutan",
                    "non_industrialized_human","industrialized_human")
physeq16s_Description_Mean$Sample <- factor(physeq16s_Description_Mean$Sample, levels =sample_ordered)

#FigureS10A
plot_phyla_abund_site <- physeq16s_DescriptionSite_Mean %>%
  filter(Phylum %in% Phylum_over01) %>%
  group_by(Sample,Phylum) %>% #group all ASVs within a Phylum
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=Phylum, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Set3") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance") +
  scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))

#FigureS10B
plot_gen_enriched_site <- physeq16s_DescriptionSite_Mean %>%
  filter(Genus %in% Genera_enriched) %>%
  group_by(Sample,Genus) %>% #group all ASVs within a Genus
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=Genus, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Dark2") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance") +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4))

#FigureS10C
plot_gen_depleted_site <- physeq16s_DescriptionSite_Mean %>%
  filter(Genus %in% Genera_depleted) %>%
  group_by(Sample,Genus) %>% #group all ASVs within a Genus
  summarise(Abundance = sum(Abundance)) %>%
  ggplot(aes(fill=Genus, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance") +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4))

#stitch plots show composition according to site together
plotleg_phyla_abund_site <- plot_grid(plot_phyla_abund_site +
                                        theme(legend.position = "none") +
                                        theme(axis.title.x=element_blank()) +
                                        theme(axis.text.x = element_blank()),get_legend(plot_phyla_abund_site),rel_widths = c(1,.4))
plotleg_gen_enriched_site <- plot_grid(plot_gen_enriched_site  +
                                         theme(legend.position = "none") +
                                         theme(axis.title.x=element_blank()) +
                                         theme(axis.text.x = element_blank()),get_legend(plot_gen_enriched_site),rel_widths = c(1,.4))
plotleg_gen_depleted_site <- plot_grid(plot_gen_depleted_site +
                                         theme(legend.position = "none"),get_legend(plot_gen_depleted_site),rel_widths = c(1,.4))
composition_barplot_site <- plot_grid(plotleg_phyla_abund_site,plotleg_gen_enriched_site,plotleg_gen_depleted_site,
                                      ncol=1,
                                      rel_heights = c(1,1,2.25),
                                      labels=c('A','B','C'))
ggsave(composition_barplot_site,file=file.path(figure_outdir,'FigS10_composition_barplot_site.pdf'),height=10,width=10)

#ALPHA DIVERSITY
calc_alpha_div <- function(physeq){
  alpha = estimate_richness(physeq, measures =  c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #calc standard alpha div metrics in phyloseq package
  pd = picante::pd(samp = t(otu_table(physeq)), tree = phy_tree(physeq), include.root = F) #calculate Faith's phylogenetic diversity
  metadata = data.frame(sample_data(physeq)) %>% select(dataset,Description_country_zoo, #extract subset of metadata
                                                        captivity_status,Description)
  alphadiv <- cbind(metadata,alpha,pd)  #merge metadata and alpha diversity estimates
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

plot_alpha_figure <- function(alpha_div_table, category, metric, color_vec) {
  #plot alpha diversity with formatting
  alpha_plot <- ggplot(alpha_div_table, aes_string(x = category, y = metric, fill=category))  +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "none") +
    ylab(paste0(metric," Diversity")) +
    scale_fill_manual(values = color_vec)
  return(alpha_plot)
}

description_site_colors <- c("indianred2","tan1","tan1","darkolivegreen3","darkolivegreen3","darkolivegreen3",
                             "plum3","plum3","plum3","plum3","blue","blue","skyblue2","skyblue2",
                             "red2","darkorange2","darkgreen")

kruskal_alphadiversity <- function(TaxName,group1,group2,df) {
  #runs kruskal-wallis and retrieves group summary stats 
  df <- df %>% filter(captivity_status == group1 | captivity_status == group2)
  group1_mean = df %>% filter(captivity_status==group1|Description==group1) %>% summarise(group1_mean = mean(Observed)) 
  group2_mean = df %>% filter(captivity_status==group2|Description==group2) %>% summarise(group2_mean = mean(Observed))
  group1_n = df %>% filter(captivity_status==group1|Description==group1) %>% nrow()
  group2_n = df %>% filter(captivity_status==group2|Description==group2) %>% nrow()
  if (group1_mean>0|group2_mean>0){ #check to make sure values in wild and captive apes
    kruskal <- df %>% 
      kruskal_test(Observed ~ captivity_status) %>% 
      select(df,statistic,p) %>% 
      unlist() } else {
        kruskal <- c(NA,NA,NA)}
  df.kruskal <- unlist(c(TaxName,group1,group2,group1_n,group2_n,group1_mean,group2_mean,kruskal))
  names(df.kruskal) <- c('TaxName','group1','group2','group1_n','group2_n',
                         'group1_mean','group2_mean','df','statistic','p')
  return(df.kruskal)
}

print('Comparing # of All ASVs in captive vs. wild apes of all host species')
alphadiv_alltaxa <- calc_alpha_div(physeq16s)
plot_alphadiv_allASVs <- plot_alpha_figure(alphadiv_alltaxa,'Description_country_zoo','Observed',description_site_colors) +
  theme(axis.text.x = element_blank()) +
  ylab('All ASVs')
(kruskal_alphadiv_allPhyla <- kruskal_alphadiversity('All Phyla','captive ape','wild ape',alphadiv_alltaxa))

#Plot alpha diversity for each of the major phyla
print('Comparing # of Actinobacteria ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Actinobacteria')
alphadiv_Actinobacteria <- calc_alpha_div(phylum_physeq16s)
plot_alphadiv_Actinobacteria <- plot_alpha_figure(alphadiv_Actinobacteria,'Description_country_zoo','Observed',description_site_colors)+
  theme(axis.text.x = element_blank()) +
  ylab('Actinobacteria ASVs')
(kruskal_alphadiv_Actinobacteria <- kruskal_alphadiversity('Actinobacteria','captive ape','wild ape',alphadiv_Actinobacteria))

print('Comparing # of Bacteroidetes ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Bacteroidetes')
alphadiv_Bacteroidetes <- calc_alpha_div(phylum_physeq16s)
plot_alphadiv_Bacteroidetes <- plot_alpha_figure(alphadiv_Bacteroidetes,'Description_country_zoo','Observed',description_site_colors)+
  theme(axis.text.x = element_blank()) +
  ylab('Bacteroidetes ASVs')
(kruskal_alphadiv_Bacteroidetes <- kruskal_alphadiversity('Bacteroidetes','captive ape','wild ape',alphadiv_Bacteroidetes))

print('Comparing # of Firmicutes ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Firmicutes')
alphadiv_Firmicutes <- calc_alpha_div(phylum_physeq16s)
plot_alphadiv_Firmicutes <- plot_alpha_figure(alphadiv_Firmicutes,'Description_country_zoo','Observed',description_site_colors)+
  theme(axis.text.x = element_blank()) +
  ylab('Firmicutes ASVs')
(kruskal_alphadiv_Firmicutes <- kruskal_alphadiversity('Firmicutes','captive ape','wild ape',alphadiv_Firmicutes))

print('Comparing # of Proteobacteria ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Proteobacteria')
alphadiv_Proteobacteria <- calc_alpha_div(phylum_physeq16s)
plot_alphadiv_Proteobacteria <- plot_alpha_figure(alphadiv_Proteobacteria,'Description_country_zoo','Observed',description_site_colors)+
  ylab('Proteobacteria ASVs')
(kruskal_alphadiv_Proteobacteria <- kruskal_alphadiversity('Proteobacteria','captive ape','wild ape',alphadiv_Proteobacteria))

#stitch together all alpha diversity plots
(alpha_all_phlya <- plot_grid(plot_alphadiv_allASVs,
                             plot_alphadiv_Actinobacteria,
                             plot_alphadiv_Bacteroidetes,
                             plot_alphadiv_Firmicutes,
                             plot_alphadiv_Proteobacteria,
                             ncol=1, rel_heights = c(1,1,1,1,2)))
ggsave(alpha_all_phlya,file=file.path(figure_outdir,'FigS7_alpha_all_phlya.pdf'),height=15,width=6)

#merge kruskal-wallis test results and apply pvalue correction
kruskal_alphadiv <- data.frame(rbind(kruskal_alphadiv_allPhyla,
      kruskal_alphadiv_Actinobacteria,
      kruskal_alphadiv_Bacteroidetes,
      kruskal_alphadiv_Firmicutes,
      kruskal_alphadiv_Proteobacteria))
kruskal_alphadiv <- kruskal_alphadiv %>% 
                      mutate(p = as.numeric(p),
                             statistic=round(as.numeric(statistic),1),
                             group1_mean=round(as.numeric(group1_mean),4),
                             group2_mean=round(as.numeric(group2_mean),4))
kruskal_alphadiv$p_adj_bonferroni <- round(p.adjust(kruskal_alphadiv$p,method = "bonferroni"),4)
kruskal_alphadiv$p <- round(kruskal_alphadiv$p,4)
write.table(as.data.frame(kruskal_alphadiv),file=file.path(table_outdir,'Table_kruskal_wallis_alphadiv_FigureS7.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#OCCUPANCY/RELATIVE ABUNDANCE
asv <- otu_table(physeq16s)
asv <- asv/colSums(asv) #relative abundance
asv_PA <- 1*((asv>0)==1)    # presence-absence data
asv_occ <- rowSums(asv_PA)  # occupancy calculation
asv_rel <- apply(decostand(asv, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- add_rownames(as.data.frame(cbind(asv_occ, asv_rel)),'ASV') # combining occupancy and abundance data frame

#visualize occupancy and relative abundance based on whether ASV is host-restricted or mixed-host
occ_abund_HRcat_filt <- occ_abun %>% 
  left_join(HR_table,on='ASV') %>% #add HRtype info
  filter(asv_occ < 200) #remove high occupancy mixed-host ASVs
(plot_occ_abund_HRcat <- ggplot(occ_abund_HRcat_filt, aes(x =log(asv_rel), y = asv_occ,color=HR_cat)) +
    geom_point() +
    geom_smooth(method="gam", aes(fill=HR_cat),formula = y ~ s(x, bs = "cs")) +
    theme_cowplot() +
    ylab('samples observed')  +
    xlab('log mean relative abundance'))
ggsave(plot_occ_abund_HRcat, file = file.path(figure_outdir,'FigS3_occ_abund_HRcat.pdf'))

print('log mean rel abund vs. samples observed model fit')
(FigureS3_models <- occ_abund_HRcat_filt %>%
  select(HR_cat,asv_occ,asv_rel) %>% nest(-HR_cat) %>%
  mutate(fit = purrr::map(data, ~mgcv::gam(asv_occ ~ s(log(asv_rel), bs = "cs"), data = .)),
         results = purrr::map(fit, glance),
         R.square = purrr::map_dbl(fit, ~ summary(.)$r.sq)) %>%
  unnest(results) %>%
  select(-data, -fit))
write.table(as.data.frame(FigureS3_models),file=file.path(table_outdir,'Table_model_summary_FigureS3.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#TAXONOMY OF HOST-RESTRICTED ASV
print('distribution of ASVs by HR category and taxonomy')
physeq16s_HR_tax <- data.frame(tax_table(physeq16s)) %>%
  group_by(Phylum,HR_cat) %>%
  tally()  %>% #summarize number of ASVs by HR category and phylum
  filter(Phylum %in% Phylum_over01) #subset to phylum reaching >1% in a sample type
(plot_physeq16s_HR_tax <- ggplot(physeq16s_HR_tax, aes(fill=Phylum, y=n, x=HR_cat)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette="Set3") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(legend.position = "none") +
  ylab("Proportion of ASVs"))
ggsave(plot_physeq16s_HR_tax,
       file = file.path(figure_outdir,'FigS4_HRtype_ASV_Taxonomy.pdf'),
       width=4)
#chi-squared
physeq16s_HR_tax_chi <- physeq16s_HR_tax %>%
  spread(key=HR_cat,value=n,fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames(var='Phylum')
(chisq <- chisq.test(physeq16s_HR_tax_chi))
chisq.res <- c(chisq$statistic,chisq$parameter,chisq$p.value)
names(chisq.res) <- c('statistic','df','p.value')
write.table(chisq.res,file=file.path(table_outdir,'Chisq_FigureS4.txt'),sep='\t',col.names=F,quote=F,row.names = TRUE)

print('relative contributions of phyla to statistic')
contrib <- 100*chisq$residuals^2/chisq$statistic
contrib <-round(contrib, 3)
write.table(contrib,file=file.path(table_outdir,'Chisq_contributions_FigureS4.txt'),sep='\t',quote=F,row.names = TRUE,col.names=NA)

#distribution of captive-ape ASVs across enclosures/host species/sites
CP_ASVs<- HR_table %>% filter(CP_pres=='True') #subset to only ASVs present in captive apes
(plot_CP_ASVs <- ggplot(CP_ASVs, aes(numEnclosure, fill=multi_site_sp)) +
    geom_bar() +
    ylab('# of 16S ASVs')+
    xlab('# of enclosures')+
    theme_bw()+
    scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a'))+
    scale_x_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10))+
    facet_wrap(~HR_cat))
ggsave(plot_CP_ASVs,file=file.path(figure_outdir,'FigS5_captive_ASVs_distribution_enclosures.pdf'),height=4)

#PERCENT SHARED ASVS
dist <- phyloseq::distance(physeq16s, method = "jaccard", binary=TRUE)
distM <- as.matrix(dist)
df <- data.frame(t(combn(colnames(distM),2)), dist=t(distM)[lower.tri(distM)]) #extract pairwise comparisons
colnames(df) <- c('ind1','ind2','dist')
site_species <- data.frame(metadata) %>% #pare down metadata
  select(captivity_status,site,Description) %>% 
  rownames_to_column(var='SampleID')
df <- df %>% 
        left_join(site_species,by=c('ind1'='SampleID')) %>%  #add metadata for individual1 and individual2
        left_join(site_species,by=c('ind2'='SampleID'),suffix = c(".ind1", ".ind2"))
df <- df %>%
  mutate(site_comp = if_else(site.ind1 == site.ind2,'same_site','diff_site')) %>% #determine if from same site
  mutate(species_comp = if_else(Description.ind1 == Description.ind2,'same_species','diff_species')) %>% #determine if same host species
  mutate(Description.ind1 = sapply(str_split(Description.ind1,pattern = '_'), "[", 2)) %>% 
  mutate(Description.ind2 = sapply(str_split(Description.ind2,pattern = '_'), "[", 2)) %>%
  mutate(species_site_comp = paste(species_comp,site_comp,sep='_')) %>%
  mutate(similarity = 1-dist) %>%
  filter(captivity_status.ind1 == 'captive' & captivity_status.ind2 == 'captive')

(Fig4_dist_site_sp_comp <- ggplot(df, aes(x=species_site_comp, y=similarity,fill=species_site_comp)) +
    geom_violin()+
    theme_bw()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared ASVs'))
ggsave(Fig4_dist_site_sp_comp,file=file.path(figure_outdir,'Fig4_dist_site_sp_comp.pdf'))

print('running stats for figure 4, comparing sorenson similarity indexes among pairwise comparisons')
perm_ttest<-function(group1,group2,df){
  group1_values <- df$dist[df$species_site_comp==group1]
  group2_values <- df$dist[df$species_site_comp==group2]
  ttest <- perm.t.test(group1_values,group2_values)
  res <- c(group1,group2,length(group1_values),length(group2_values),
           ttest$estimate[1],abs(ttest$estimate[2]),ttest$statistic,ttest$p.value)
  names(res) <- c('group1','group2','group1_n','group2_n',
                  'group1_mean','group2_mean','statistic','p.value')
  return(res)
}

Fig4_stats_table <- data.frame(rbind(
  perm_ttest("diff_species_diff_site","diff_species_same_site",df),
  perm_ttest("diff_species_diff_site","same_species_diff_site",df),
  perm_ttest("diff_species_diff_site","same_species_same_site",df),
  perm_ttest("diff_species_same_site","same_species_diff_site",df),
  perm_ttest("diff_species_same_site","same_species_same_site",df),
  perm_ttest("same_species_diff_site","same_species_same_site",df)))

Fig4_stats_table <- Fig4_stats_table %>%
                          mutate(p.value = as.numeric(p.value),
                                 statistic=round(as.numeric(statistic),1),
                                 group1_mean=round(as.numeric(group1_mean),4),
                                 group2_mean=round(as.numeric(group2_mean),4))
Fig4_stats_table$p_adj_bonferroni <- round(p.adjust(Fig3_stats_table$p.value,method = "bonferroni"),4)
Fig4_stats_table$p.value <- round(Fig3_stats_table$p.value,4)
write.table(as.data.frame(Fig4_stats_table),file=file.path(table_outdir,'Table_permttest_Figure4.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#Figure S6
df_same_site <- df %>% filter(species_site_comp=='diff_species_same_site')
df_same_site <- df_same_site %>%
  mutate(within_sites_comp=paste0(site.ind1,'_',Description.ind1,'_vs_',Description.ind2)) %>%
  mutate(within_sites_comp = gsub("como_zoo_orangutan_vs_gorilla","como_zoo_gorilla_vs_orangutan", within_sites_comp))

plot_shared_ASVs_wn_site <- ggplot(df_same_site, aes(x=within_sites_comp, y=similarity,fill=within_sites_comp)) +
  geom_violin()+
  scale_fill_manual(values=c('orange2','orange2','orange2','green3',
                             'skyblue','skyblue','skyblue'))+
  theme_bw()+
  ylab('Proportion of shared ASVs')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot_shared_ASVs_wn_site,
       file = file.path(figure_outdir,'FigS6_shared_ASVs_wn_sites.pdf'))

#exclude chimpvsgorilla houston from the analysis
#to see how much these samples contribute to the analysis
diff_species_diff_site <- df$dist[df$species_site_comp=='diff_species_diff_site']
diff_species_same_site_excludeHOUS <- df_same_site$dist[df_same_site$within_sites_comp!='houston_zoo_chimp_vs_gorilla']
perm.t.test(diff_species_diff_site,diff_species_same_site_excludeHOUS)
#still significant, but t-statistic drops from 16.3 to 6.9
#when excluding Houston zoo comparisons
#system('rm Rplots.pdf') #remove plots generated by capture-output cmd 
