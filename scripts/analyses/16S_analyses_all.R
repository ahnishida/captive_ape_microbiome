library("phyloseq")
library("tidyverse")
library("picante")
library("PMCMRplus")
library("cowplot")
library("reshape2")
library("RVAideMemoire")
library("RColorBrewer")
library('rstatix')

setwd('/Volumes/AHN/captive_ape_microbiome/') #SET WORKING DIRECTORY
#INPUTS
physeq16s <- readRDS('results/16s/inputs/phyloseq_rare10000.rds') #read in physeq
HR_table <- read.table('results/16s/analyses/tables/16S_ASVs_summary.txt',sep='\t',header=T)
metadata <- sample_data(physeq16s) #format metadata

#OUTDIR
table_outdir <- 'results/16s/analyses/tables' #specify output folder for tables
dir.create(table_outdir,recursive=TRUE)
figure_outdir <-'results/16s/analyses/figures' #specify output folder for figures
dir.create(figure_outdir,recursive=TRUE)

#set colors for groups
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
capture.output(plot_nmds_bray <- ordinate_NMDS(physeq16s,physeq16s_bray),file='nul') #Figure1
ggsave(plot_nmds_bray, file = file.path(figure_outdir,'Fig1_nMDS_bray.pdf'),width=6,height=5)
capture.output(plot_nmds_jaccard <- ordinate_NMDS(physeq16s,physeq16s_jaccard),file='nul')
capture.output(plot_nmds_wunifrac <- ordinate_NMDS(physeq16s,physeq16s_wunifrac),file='nul')
capture.output(plot_nmds_uunifrac <- ordinate_NMDS(physeq16s,physeq16s_uunifrac),file='nul')
plot_nmds_othermetrics <- plot_grid(plot_nmds_jaccard,
                                    plot_nmds_wunifrac,
                                    plot_nmds_uunifrac,
                                    ncol=1,
                                    labels = "AUTO")
ggsave(plot_nmds_othermetrics, file = file.path(figure_outdir,'FigS2_nMDS_othermetrics.pdf'),width=6,height=14)

#BETADISPER and PERMANOVA - Description all groups
run_betadisper_permanova_dist <- function(physeq,physeq_dist){
  metadata = as.data.frame(as.matrix(sample_data(physeq))) %>% #extract metadata
    rownames_to_column(var="X.SampleID")
  print('testing homogeneity of sample groups')
  beta <- betadisper(physeq_dist, metadata$Description) #run betadisper
  beta_tab <- permutest(beta)$tab #extract table
  print(beta_tab)
  beta_p <- beta_tab$P[[1]]
  print('permanova')
  perm <- adonis(physeq_dist ~ Description, data = metadata) #run permanova
  perm_tab <- as.data.frame(perm$aov.tab) #extract table
  print(perm_tab)
  perm_pseudof <- perm_tab$F.Model[[1]]
  perm_R2 <- perm_tab$R2[[1]]
  perm_pvalue <- perm_tab$Pr[[1]]
  beta_perm_res <- c(beta_p,perm_pseudof,perm_R2,perm_pvalue) #combine betadisper and permanova
  return(beta_perm_res)
}

print('BETADISPER and PERMANOVA')
print('bray-curtis')
(beta_perm_bray <- run_betadisper_permanova_dist(physeq16s,physeq16s_bray))
print('jaccard')
(beta_perm_jaccard <- run_betadisper_permanova_dist(physeq16s,physeq16s_jaccard))
print('weighted unifrac')
(beta_perm_wunifrac <- run_betadisper_permanova_dist(physeq16s,physeq16s_wunifrac))
print('unweighted unifrac')
(beta_perm_uunifrac <- run_betadisper_permanova_dist(physeq16s,physeq16s_uunifrac))

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
  pw <- pw_comps[i,]
  print(pw) #pairwise comparison
  physeq_pw <- subset_samples(physeq16s,Description==pw$group1|Description==pw$group2) #subset physeq to two groups
  physeq_pw_dist <- phyloseq::distance(physeq_pw,'bray') #calculate distance
  capture.output(beta_perm_res <- run_betadisper_permanova_dist(physeq_pw,physeq_pw_dist),file='nul') #run betadisper and permanova
  res <- c(pw$group1,pw$group2,beta_perm_res) #generate new line with group names and results
  pw_df <- rbind(pw_df,res) #bind new row to dataframe
  }
colnames(pw_df) <- c('group1','group2','betadisper_Pr(>F)','PERMANOVA_F.Model','PERMANOVA_R2','PERMANOVA_Pr(>F)')
write.table(pw_df,file=file.path(table_outdir,'Table_pw_betadisper_permanova_bray.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#ADD HR INFO to PHYSEQ TAX_TABLE
new_tax_table <- data.frame(tax_table(physeq16s)) #add HRcat and type to taxonomy
new_tax_table$HR_cat <- HR_table$HR_cat
new_tax_table$HR_type <- HR_table$HR_type
tax_table(physeq16s) <- as.matrix(new_tax_table)

#Make new column in metadata that is a combo of description/zoo site
metadata <- sample_data(physeq16s) #format metadata with new column
metadata$Description_country_zoo <- paste(metadata$Description,metadata$country_zoo,sep='_')
sample_data(physeq16s) <- metadata

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
Phylum_over01 <- physeq16s_Description_Mean %>% 
  group_by(Phylum) %>%
  summarise(max_abundance=max(Abundance)) %>%
  as.data.frame() %>%
  filter(max_abundance>.01)

print('summary of the mean relative abundance of bacterial phylum across sample groups')
(physeq16s_Description_Mean %>% 
    filter(Phylum %in% Phylum_over01$Phylum) %>%
    group_by(Phylum,Sample) %>% 
    summarise(mean_abundance=sum(Abundance)) %>%
    as.data.frame() %>% 
    spread(key = Phylum, value=mean_abundance))
print('summary of the mean relative abundance of HR ASVs across sample groups')
(physeq16s_Description_Mean %>% 
    group_by(HR_cat,Sample) %>% 
    summarise(mean_abundance=sum(Abundance)) %>%
    as.data.frame() %>% 
    spread(key = HR_cat, value=mean_abundance))
(physeq16s_Description_Mean %>% 
  group_by(HR_type,Sample) %>% 
  summarise(mean_abundance=sum(Abundance)) %>%
  as.data.frame() %>% 
  spread(key = HR_type, value=mean_abundance))


plot_phyla_abund <- physeq16s_Description_Mean %>%
  group_by(Sample,Phylum) %>% #group all ASVs within a HR_type
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  filter(Phylum %in% Phylum_over01$Phylum) %>%
  ggplot(aes(fill=Phylum, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Set3") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))

#COMPOSITION BY PHYLA - STATS
physeq16s_stats <- physeq16s #make copy of physeq
otu_table(physeq16s_stats) <- otu_table(physeq16s_stats)/sample_sums(physeq16s_stats) #format count to relative abundance
physeq16s_melt <- psmelt(physeq16s_stats) #combines taxa table, metadata, and otu table into df

kw_abund_Taxa <- function(TaxLevel,TaxName,nComp){
  #runs kruskal-wallis and posthoc dunn test on abundance of bacterial phylum across group descriptions
  physeq16s_Taxa <- physeq16s_melt %>%
    filter(!!sym(TaxLevel)==TaxName) %>% #subset to the just taxa of interest
    group_by(Sample,captivity_status,Description) %>%  #group all ASVs within an individual, keep group metadata
    summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
    as.data.frame()
   #get captive vs wild  means
  cpwd_means <- physeq16s_Taxa %>%
    filter(captivity_status=='wild'|captivity_status=='captive') %>%
    group_by(captivity_status) %>%
    summarise(Mean = mean(Abundance)) %>%
    as.data.frame() %>%
    column_to_rownames(var='captivity_status') %>%
    t()

  #get description means
  group_means <- physeq16s_Taxa %>%
    group_by(Description) %>%
    summarise(Mean = mean(Abundance)) %>%
    as.data.frame() %>%
    column_to_rownames(var='Description') %>%
    t()

  #captive vs wild all species
  cp_wd <- physeq16s_Taxa %>% filter(captivity_status=='captive'|captivity_status=='wild')
  if (max(cp_wd$Abundance)>0){ #check to make sure values in wild and captive apes
    res.kruskal <- cp_wd  %>% kruskal_test(Abundance ~ captivity_status)
    (cp_wd_kw_adj <- p.adjust(res.kruskal$p,method = "fdr",n = nComp))} else {cp_wd_kw_adj <- NA}
  
  #captive vs wild bonobo
  bonobo_comp <- physeq16s_Taxa %>% filter(Description=='captive_bonobo'|Description=='wild_bonobo')
  if (max(bonobo_comp$Abundance)>0){
    res.kruskal <- bonobo_comp  %>% kruskal_test(Abundance ~ captivity_status)
    (bonobo_kw_pval_adj <- p.adjust(res.kruskal$p,method = "fdr",n = nComp))} else {bonobo_kw_pval_adj <- NA}
  
  #captive vs wild chimp
  chimp_comp <- physeq16s_Taxa %>% filter(Description=='captive_chimp'|Description=='wild_chimp')
  if (max(chimp_comp$Abundance)>0){
    res.kruskal <- chimp_comp %>% kruskal_test(Abundance ~ captivity_status)
    (chimp_kw_pval_adj <- p.adjust(res.kruskal$p,method = "fdr",n = nComp))} else {chimp_kw_pval_adj <- NA}
  
  #captive vs wild gorilla
  gorilla_comp <- physeq16s_Taxa %>% filter(Description=='captive_gorilla'|Description=='wild_gorilla')
  if (max(gorilla_comp$Abundance)>0){
    res.kruskal <- gorilla_comp %>% kruskal_test(Abundance ~ captivity_status)
    (gorilla_kw_pval_adj <- p.adjust(res.kruskal$p,method = "fdr",n = nComp))} else {gorilla_kw_pval_adj <- NA}
  
  res <- c(TaxName,cpwd_means,group_means,cp_wd_kw_adj,bonobo_kw_pval_adj,chimp_kw_pval_adj,gorilla_kw_pval_adj)
  names(res) <- c(TaxLevel,paste0(colnames(cpwd_means),'_mean'),paste0(colnames(group_means),'_mean'),'captive_vs_wild_allspecies','captive_vs_wild_bonobo','captive_vs_wild_chimp','captive_vs_wild_gorilla')
  return(res)
}

print('kruskal-wallis differential abundance test')
kw_abund_Taxa('Phylum','Actinobacteria',4) #test on a single phylum
kw_abund_Taxa('Genus','Bifidobacterium',4) #test on a single genus

#run Kruskal-Wallis and post hoc on all phyla over 1%
Phyla_nComp <- length(Phylum_over01$Phylum)*4 #function run 4 comparisons for each taxa
kw_abund_Phyla <- function(TaxName){
  kw_abund_Taxa('Phylum',TaxName,Phyla_nComp)
}
kw_phylum <- lapply(Phylum_over01$Phylum,kw_abund_Phyla)
kw_phylum_df <- as.data.frame(do.call(rbind,kw_phylum)) #combine into dataframe
num_col <- colnames(kw_phylum_df)[2:length(colnames(kw_phylum_df))] #numeric columns
kw_phylum_df <- kw_phylum_df %>% mutate_each_(funs(as.numeric), num_col) #change from character to numeric
write.table(kw_phylum_df,file=file.path(table_outdir,'Table_kruskal_wallis_phyla.txt'),sep='\t',row.names = F,col.names=T,quote=F)

print(kw_phylum_df %>% select(Phylum,captive_mean,wild_mean,
                              captive_vs_wild_allspecies,captive_vs_wild_bonobo,
                              captive_vs_wild_chimp,captive_vs_wild_gorilla))
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
plot_phyla_abund_HR <- physeq16s_Description_Mean %>%
  group_by(Sample,HR_type) %>% #group all ASVs within a HR_type
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=HR_type, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = HRpalette) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))


#COMPOSITION BY HR_TYPE - STATS
print('running statistics compare relative abundance of host-restricted between sample groups')

kw_dunn <- function(df,nComp) {
  #runs Kruskal Wallis and if significant runs dunns posthoc
  res.kruskal <- df %>% kruskal_test(Abundance ~ captivity_status)
  #print(res.kruskal)
  pval <- p.adjust(res.kruskal$p,method = "fdr",n = nComp)
  print(c(pval,'pval_adjusted'))
  if (length(unique(df$captivity_status))>2 & pval<.05){
    df$captivity_status <- as.factor(df$captivity_status)
    print(kwAllPairsDunnTest(Abundance ~ captivity_status, df,
                             p.adjust.method = "fdr")) }
}

HR_type_stats <-function(physeq_melt){
  
  df <- physeq_melt %>%
    group_by(Sample,captivity_status,Description,HR_type) %>%  #group all ASVs within an individual, keep group metadata
    summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
    as.data.frame()
  df$captivity_status <- recode(df$captivity_status,non_industrialized_human='human') 
  df$captivity_status <- recode(df$captivity_status,industrialized_human='human')
  
  print('comparing relative abundance of HR_human ASVs, humans vs. all captive ape species')
  HR_human <- df %>% #remove wild apes
    filter(HR_type == 'HR_human' & captivity_status != 'wild')
  #recode non-industrialized and industrialized humans as 
  kw_dunn(HR_human,6)

  print('comparing relative abundance of HR_wild_gorilla ASVs, wild gorillas vs. all captive ape species')
  HR_wild_gorilla  <- df %>%
    filter(HR_type == 'HR_wild_gorilla') %>%
    filter(captivity_status == "captive" | Description == "wild_gorilla")
  kw_dunn(HR_wild_gorilla,6)

  print('comparing relative abundance of HR_wild_chimp ASVs, wild chimps vs. all captive ape species')
  HR_wild_chimp  <- df %>%
    filter(HR_type == 'HR_wild_chimp') %>%
    filter(captivity_status == "captive" | Description == "wild_chimp")
  kw_dunn(HR_wild_chimp,6)

  print('comparing relative abundance of HR_wild_bonobo ASVs, wild bonobos vs. all captive ape species')
  HR_wild_bonobo  <- df %>%
    filter(HR_type == 'HR_wild_bonobo') %>%
    filter(captivity_status == "captive"| Description == "wild_bonobo")
  kw_dunn(HR_wild_chimp,6)

  print('comparing relative abundance of MX_wild_apes ASVs, wild apes vs. captive apes')
  MX_wild_apes  <- df %>%
    filter(HR_type == 'MX_wild_apes' & captivity_status != 'human') 
  kw_dunn(MX_wild_apes,6)

  print('comparing relative abundance of MX_human_wild_apes ASVs, wild apes vs. captive apes. vs humans')
  MX_human_wild_apes  <- df %>%
    filter(HR_type == 'MX_human_wild_apes')
  kw_dunn(MX_human_wild_apes,6)

}
HR_type_stats(physeq16s_melt)

#BACTERIAL GENERA DIFF ABUNDANCE IN WILD AND CAPTIVE APES
Genera <- unique(physeq16s_melt$Genus) #list of unique genera
Genera <- Genera[!is.na(Genera)] #remove NA from list
Genera_nComp <- length(Genera)*4
kw_abund_Genera <- function(TaxName){
  kw_abund_Taxa('Genus',TaxName,Genera_nComp)
}

kw_Genera_df <- lapply(Genera,kw_abund_Genera)
kw_Genera_df <- as.data.frame(do.call(rbind, kw_Genera_df)) #combine into dataframe
num_col <- colnames(kw_Genera_df)[2:length(colnames(kw_Genera_df))] #numeric columns
kw_Genera_df <- kw_Genera_df %>% mutate_each_(funs(as.numeric), num_col) #change from character to numeric
write.table(kw_Genera_df,file=file.path(table_outdir,'Table_kruskal_wallis_genera.txt'),sep='\t',row.names = F,col.names=T,quote=F)

#Identify bacterial genera that show parallel enrichment in all captive ape species
Genera_enriched_captive <- kw_Genera_df %>% filter(captive_vs_wild_allspecies < .05) %>%
                                  filter(captive_mean > wild_mean) %>%
                                  filter(captive_vs_wild_bonobo < .05) %>%
                                  filter(captive_vs_wild_chimp < .05) %>%
                                  filter(captive_vs_wild_gorilla < .05) %>%
                                  filter(captive_bonobo_mean > wild_bonobo_mean) %>%
                                  filter(captive_chimp_mean > wild_chimp_mean) %>%
                                  filter(captive_gorilla_mean > wild_gorilla_mean)
Genera_enriched_captive_over01 <- Genera_enriched_captive %>%
  filter(captive_bonobo_mean>.01&captive_chimp_mean>.01&captive_gorilla_mean>.01)
print('Bacterial genera enriched in all captive ape species present at >1% abundance')
print(Genera_enriched_captive_over01 %>% select(Genus,captive_mean,wild_mean,
                                                    captive_vs_wild_allspecies,captive_vs_wild_bonobo,
                                                    captive_vs_wild_chimp,captive_vs_wild_gorilla))
#Identify bacterial genera depleted in captive vs. wild
Genera_depleted_captive_onecomp <- kw_Genera_df %>% 
  filter(captive_mean < wild_mean & captive_vs_wild_allspecies < .05) %>%
  filter(
  (as.numeric(captive_vs_wild_bonobo) < .05 & captive_bonobo_mean < wild_bonobo_mean & wild_bonobo_mean>.02) |
  (captive_vs_wild_chimp < .05 & captive_chimp_mean < wild_chimp_mean & wild_chimp_mean> .02) |
  (captive_vs_wild_gorilla < .05 & captive_gorilla_mean < wild_gorilla_mean & wild_gorilla_mean > .02))
print('Bacterial genera depleted in captive apes vs wild apes in one species comparison, present at >2% abundance')
print(Genera_depleted_captive_onecomp %>% select(Genus,captive_mean,wild_mean,
                                                captive_vs_wild_allspecies,captive_vs_wild_bonobo,
                                                captive_vs_wild_chimp,captive_vs_wild_gorilla))

print('run HR stats on genera enriched/depleted')
print('shows same trends across all ASVs')
print('HR stats for genera enriched in captive apes')
physeq16s_melt_Genera_enriched <- physeq16s_melt %>%
  filter(Genus %in% Genera_enriched_captive_over01$Genus)
HR_type_stats(physeq16s_melt_Genera_enriched)
print('HR stats for genera depleted in captive apes')
physeq16s_melt_Genera_depleted <- physeq16s_melt %>%
  filter(Genus %in% Genera_depleted_captive_onecomp$Genus)
HR_type_stats(physeq16s_melt_Genera_depleted)

#Visualize bacterial genera enriched in captivity
(plot_genus_enriched_abund <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_enriched_captive_over01$Genus) %>%
  group_by(Sample,Genus) %>% #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=Genus, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Dark2") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

(plot_genus_enriched_abund_HRtype <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_enriched_captive_over01$Genus) %>%
  group_by(Sample,HR_type) %>%  #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=HR_type, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = HRpalette) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance") +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

#Visualize bacterial genera depleted in captivity
colors = c(brewer.pal(name="Paired", n = 12),brewer.pal(name="Dark2", n = 8))
(plot_genus_deplete_abund <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_depleted_captive_onecomp$Genus) %>%
  group_by(Sample,Genus) %>% #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=Genus, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

(plot_genus_deplete_abund_HRtype <- physeq16s_Description_Mean %>%
  filter(Genus %in% Genera_depleted_captive_onecomp$Genus) %>%
  group_by(Sample,HR_type) %>% #group all ASVs within an individual, keep group metadata
  summarise(Abundance = sum(Abundance)) %>% #merge ASVs to get their sum
  ggplot(aes(fill=HR_type, y=Abundance, x=Sample)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = HRpalette) +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90,vjust =.4,hjust=1)) +
  ylab("Relative abundance")  +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4)))

#Figure 2
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
ggsave(composition_taxa_HR,file=file.path(figure_outdir,'Fig2_composition_taxa_HR.pdf'),height=11,width=12)

#Figure S3
plot_genus_deplete_abund_combined <- plot_grid(plot_genus_deplete_abund + theme(legend.position = "none"),get_legend(plot_genus_deplete_abund),
                                               plot_genus_deplete_abund_HRtype + theme(legend.position = "none"),get_legend(plot_genus_deplete_abund_HRtype),
                                               rel_widths = c(1,.5,1,.5))
ggsave(plot_genus_deplete_abund_combined,file=file.path(figure_outdir,'FigS3_genera_depleted_captive.pdf'),height=13,width=8)


#COMPOSITION BARPLOT BY SITE
physeq16s_DescriptionSite_Mean <-  merge_samples(physeq16s, "Description_country_zoo") #collaspe samples into groups based on species and captivity status
otu_table(physeq16s_DescriptionSite_Mean) <- otu_table(physeq16s_DescriptionSite_Mean)/sample_sums(physeq16s_DescriptionSite_Mean) #format count to relative abundance
physeq16s_DescriptionSite_Mean <- psmelt(physeq16s_DescriptionSite_Mean)
sample_ordered <- c("wild_bonobo","wild_chimp","wild_gorilla","captive_bonobo","captive_chimp",
                    "captive_gorilla","captive_orangutan",
                    "non_industrialized_human","industrialized_human")
physeq16s_Description_Mean$Sample <- factor(physeq16s_Description_Mean$Sample, levels =sample_ordered)

plot_phyla_abund_site <- physeq16s_DescriptionSite_Mean %>%
  filter(Phylum %in% Phylum_over01$Phylum) %>%
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

plot_gen_enriched_site <- physeq16s_DescriptionSite_Mean %>%
  filter(Genus %in% Genera_enriched_captive_over01$Genus) %>%
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

plot_gen_depleted_site <- physeq16s_DescriptionSite_Mean %>%
  filter(Genus %in% Genera_depleted_captive_onecomp$Genus) %>%
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

#stitch plot show composition according to site together
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
ggsave(composition_barplot_site,file=file.path(figure_outdir,'FigS4_composition_barplot_site.pdf'),height=10,width=10)

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

print('run alpha diversity statistics')
kw_dunn_alpha <- function(df,nComp) {
  #runs Kruskal Wallis and if significant runs dunns posthoc
  df <- df %>% filter(captivity_status == 'wild ape' | captivity_status == 'captive ape' )
  print(df %>% group_by(captivity_status) %>% summarise(mean(Observed)))
  (kw.res <- df %>% kruskal_test(Observed ~ captivity_status))
  pval <- p.adjust(kw.res$p,method = "fdr",n = nComp)
  print(c(pval,'pval_adjusted'))
  
 
}

description_site_colors <- c("indianred2","tan1","tan1","darkolivegreen3","darkolivegreen3","darkolivegreen3",
                             "plum3","plum3","plum3","plum3","blue","blue","skyblue2","skyblue2",
                             "red2","darkorange2","darkgreen")

alphadiv_alltaxa <- calc_alpha_div(physeq16s)
write.table(alphadiv_alltaxa,file=file.path(table_outdir,'alpha_div_16S.txt'),
            sep='\t', quote=F, row.names = F)

print('Comparing # of All ASVs in captive vs. wild apes of all host species')
plot_alphadiv_allASVs <- plot_alpha_figure(alphadiv_alltaxa,'Description_country_zoo','Observed',description_site_colors) +
  theme(axis.text.x = element_blank()) +
  ylab('All ASVs')
kw_dunn_alpha(alphadiv_alltaxa,5)

#Plot alpha diversity for each of the major phyla
print('Comparing # of Actinobacteria ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Actinobacteria')
Actinobacteria_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
alphadiv_Actinobacteria <- calc_alpha_div(Actinobacteria_physeq16s)
plot_alphadiv_Actinobacteria <- plot_alpha_figure(alphadiv_Actinobacteria,'Description_country_zoo','Observed',description_site_colors)+
  theme(axis.text.x = element_blank()) +
  ylab('Actinobacteria ASVs')
kw_dunn_alpha(alphadiv_Actinobacteria,5)

print('Comparing # of Bacteroidetes ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Bacteroidetes')
Bacteroidetes_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
alphadiv_Bacteroidetes <- calc_alpha_div(Bacteroidetes_physeq16s)
plot_alphadiv_Bacteroidetes <- plot_alpha_figure(alphadiv_Bacteroidetes,'Description_country_zoo','Observed',description_site_colors)+
  theme(axis.text.x = element_blank()) +
  ylab('Bacteroidetes ASVs')
kw_dunn_alpha(alphadiv_Bacteroidetes,5)

print('Comparing # of Firmicutes ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Firmicutes')
Firmicutes_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
alphadiv_Firmicutes <- calc_alpha_div(Firmicutes_physeq16s)
plot_alphadiv_Firmicutes <- plot_alpha_figure(alphadiv_Firmicutes,'Description_country_zoo','Observed',description_site_colors)+
  theme(axis.text.x = element_blank()) +
  ylab('Firmicutes ASVs')
kw_dunn_alpha(alphadiv_Firmicutes,5)

print('Comparing # of Proteobacteria ASVs in captive vs. wild apes of all host species')
phylum_physeq16s <- subset_taxa(physeq16s,Phylum=='Proteobacteria')
Proteobacteria_physeq16s <- prune_samples(sample_sums(phylum_physeq16s)>0, phylum_physeq16s)
alphadiv_Proteobacteria <- calc_alpha_div(Proteobacteria_physeq16s)
plot_alphadiv_Proteobacteria <- plot_alpha_figure(alphadiv_Proteobacteria,'Description_country_zoo','Observed',description_site_colors)+
  ylab('Proteobacteria ASVs')
kw_dunn_alpha(alphadiv_Proteobacteria,5)

alpha_all_phlya <- plot_grid(plot_alphadiv_allASVs,
                             plot_alphadiv_Actinobacteria,
                             plot_alphadiv_Bacteroidetes,
                             plot_alphadiv_Firmicutes,
                             plot_alphadiv_Proteobacteria,
                             ncol=1, rel_heights = c(1,1,1,1,2))
alpha_all_phlya
ggsave(alpha_all_phlya,file=file.path(figure_outdir,'FigS1_alpha_all_phlya.pdf'),height=15,width=6)

#OCCUPANCY/RELATIVE ABUNDANCE
asv <- otu_table(physeq16s)
asv <- asv/colSums(asv) #relative abundance
asv_PA <- 1*((asv>0)==1)    # presence-absence data
asv_occ <- rowSums(asv_PA)  # occupancy calculation
asv_rel <- apply(decostand(asv, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- add_rownames(as.data.frame(cbind(asv_occ, asv_rel)),'ASV') # combining occupancy and abundance data frame

#visualize occupancy and relative abundance based on whether ASV is host-restricted or mixed-host
occ_abund_HRcat <- occ_abun %>% left_join(HR_table,on='ASV')
df <- occ_abund_HRcat %>% filter('HR_cat'=='HR')
occ_abund_HRcat_filt <- occ_abund_HRcat %>% filter(asv_occ < 200)
(plot_occ_abund_HRcat <- ggplot(occ_abund_HRcat_filt, aes(x =log(asv_rel), y = asv_occ,color=HR_cat)) +
    geom_point() +
    geom_smooth(method="gam", aes(fill=HR_cat),formula = y ~ s(x, bs = "cs")) +
    theme_cowplot() +
    ylab('samples observed')  +
    xlab('log mean relative abundance'))
ggsave(plot_occ_abund_HRcat, file = file.path(figure_outdir,'FigS5_occ_abund_HRcat.pdf'))

library(broom)
print('log mean rel abund vs. samples observed model fit')
occ_abund_HRcat_filt %>%
  select(HR_cat,asv_occ,asv_rel) %>% nest(-HR_cat) %>%
  mutate(fit = map(data, ~mgcv::gam(asv_occ ~ s(log(asv_rel), bs = "cs"), data = .)),
         results = map(fit, glance),
         R.square = map_dbl(fit, ~ summary(.)$r.sq)) %>%
  unnest(results) %>%
  select(-data, -fit)

#TAXONOMY OF HOST-RESTRICTED ASV
print('distribution of ASVs by HR category and taxonomy')
table(tax_table(physeq16s)[,'HR_cat'])
data.frame(tax_table(physeq16s)) %>%
  filter(Phylum %in% Phylum_over01$Phylum) %>%
  group_by(Phylum,HR_cat) %>%
  tally() %>%
  spread(key=HR_cat,value=n)

physeq16s_HR_tax <- data.frame(tax_table(physeq16s)) %>%
  group_by(Phylum,HR_cat) %>%
  tally()  %>% #summarize number of ASVs by HR category and phylum
  filter(Phylum %in% Phylum_over01$Phylum) #subset to phylum reaching >1% in a sample type
plot_physeq16s_HR_tax <- ggplot(physeq16s_HR_tax, aes(fill=Phylum, y=n, x=HR_cat)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette="Set3") +
  theme_cowplot() +
  theme(axis.title.x=element_blank()) +
  theme(legend.position = "none") +
  ylab("Proportion of ASVs")
ggsave(plot_physeq16s_HR_tax,
       file = file.path(figure_outdir,'FigS6_HRtype_ASV_Taxonomy.pdf'),
       width=4)
#chi-squared
physeq16s_HR_tax_chi <- physeq16s_HR_tax %>%
  spread(key=HR_cat,value=n,fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames(var='Phylum')
(chisq <- chisq.test(physeq16s_HR_tax_chi))
print('relative contributions of phyla to statistic')
contrib <- 100*chisq$residuals^2/chisq$statistic
(round(contrib, 3))

#distribution of captive-ape ASVs across enclosures/host species/sites
CP_ASVs<- HR_table %>% filter(CP_pres=='True')
(plot_CP_ASVs <- ggplot(CP_ASVs, aes(numEnclosure, fill=multi_site_sp)) +
    geom_bar() +
    ylab('# of 16S ASVs')+
    xlab('# of enclosures')+
    theme_bw()+
    scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a'))+
    scale_x_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10))+
    facet_wrap(~HR_cat))
ggsave(plot_CP_ASVs,file=file.path(figure_outdir,'FigS7_captive_ASVs_distribution_enclosures.pdf'),height=4)

#PERCENT SHARED ASVS
dist <- distance(physeq16s, method = "jaccard", binary=TRUE)
distM <-as.matrix(dist)
df <- data.frame(t(combn(colnames(distM),2)), dist=t(distM)[lower.tri(distM)])
colnames(df) <- c('ind1','ind2','dist')
site_species <- data.frame(metadata) %>% select(captivity_status,site,Description) %>% rownames_to_column(var='SampleID')
df <- df %>% left_join(site_species,by=c('ind1'='SampleID')) %>% left_join(site_species,by=c('ind2'='SampleID'),suffix = c(".ind1", ".ind2"))
df <- df %>%
  mutate(site_comp = if_else(site.ind1 == site.ind2,'same_site','diff_site')) %>%
  mutate(species_comp = if_else(Description.ind1 == Description.ind2,'same_species','diff_species')) %>%
  mutate(Description.ind1 = sapply(str_split(Description.ind1,pattern = '_'), "[", 2)) %>%
  mutate(Description.ind2 = sapply(str_split(Description.ind2,pattern = '_'), "[", 2)) %>%
  mutate(species_site_comp = paste(species_comp,site_comp,sep='_')) %>%
  mutate(similarity = 1-dist) %>%
  filter(captivity_status.ind1 == 'captive' & captivity_status.ind2 == 'captive')

(Fig3_dist_site_sp_comp <- ggplot(df, aes(x=species_site_comp, y=similarity,fill=species_site_comp)) +
    geom_violin()+
    theme_bw()+
    scale_fill_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99'))+
    ylab('Proportion of shared ASVs'))
ggsave(Fig3_dist_site_sp_comp,file=file.path(figure_outdir,'Fig3_dist_site_sp_comp.pdf'))

print('running stats for figure 3, comparing sorenson similarity indexes among pairwise comparisons')
diff_species_diff_site <- df$dist[df$species_site_comp=='diff_species_diff_site']
same_species_diff_site <- df$dist[df$species_site_comp=='same_species_diff_site']
diff_species_same_site <- df$dist[df$species_site_comp=='diff_species_same_site']
same_species_same_site <- df$dist[df$species_site_comp=='same_species_same_site']

(test <- perm.t.test(diff_species_diff_site,same_species_diff_site))
print(c('fdr adj pvalue',p.adjust(test$p.value,method = "fdr",n = 6)))
(test <- perm.t.test(diff_species_diff_site,diff_species_same_site))
print(c('fdr adj pvalue',p.adjust(test$p.value,method = "fdr",n = 6)))
(test <-perm.t.test(diff_species_diff_site,same_species_same_site))
print(c('fdr adj pvalue',p.adjust(test$p.value,method = "fdr",n = 6)))
(test <-perm.t.test(diff_species_same_site,same_species_diff_site))
print(c('fdr adj pvalue',p.adjust(test$p.value,method = "fdr",n = 6)))
(test <-perm.t.test(diff_species_same_site,same_species_same_site))
print(c('fdr adj pvalue',p.adjust(test$p.value,method = "fdr",n = 6)))
(test <-perm.t.test(same_species_diff_site,same_species_same_site))
print(c('fdr adj pvalue',p.adjust(test$p.value,method = "fdr",n = 6)))

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
       file = file.path(figure_outdir,'FigS8_shared_ASVs_wn_sites.pdf'))

#exclude chimpvsgorilla houston from the analysis
#to see how much these samples contribute to the analysis
diff_species_same_site_excludeHOUS <- df_same_site$dist[df_same_site$within_sites_comp!='houston_zoo_chimp_vs_gorilla']
(test <- perm.t.test(diff_species_diff_site,diff_species_same_site_excludeHOUS))
print(c('fdr adj pvalue',p.adjust(test$p.value,method = "fdr",n = 6)))
