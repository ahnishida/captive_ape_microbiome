library("phyloseq")
library("ggplot2")
library("picante")
theme_set(theme_bw())

setwd('/Volumes/AHN/captive_ape_microbiome/results/16s')
(physeq16s <- readRDS('inputs/phyloseq_rare10000.rds'))

#ALPHA DIVERSITY
alpha = estimate_richness(physeq16s, measures = NULL)
pd = picante::pd(samp = t(otu_table(physeq16s)), tree = phy_tree(physeq16s), include.root = F)
metadata = data.frame(sample_data(physeq16s)) %>% select(dataset,genus_sp,subspecies,common_name,
                                                      site,site_code,country,
                                                      captivity_status,captivity_status_ext,Description)
alphadiv <- cbind(metadata,alpha,pd) 
alphadiv <- data.frame(alphadiv) %>% rownames_to_column(var='SampleID')
write.table(alphadiv,file='analyses/figures/Table_alpha_div_16S.txt',sep='\t',quote=F,row.names = F)

#graph FaithPD 
alphadiv$captivity_status <- as.vector(recode(alphadiv$captivity_status, wild = 'wild ape', 
                                          captive = 'captive ape',
                                          western_human = "western human",
                                          non_western_human="non western human"))
new_order = c("wild ape","captive ape","non western human","western human")
alphadiv$captivity_status <- factor(alphadiv$captivity_status, levels = new_order)
ape_palette = c("darkorange1",'deepskyblue',"goldenrod1",'mediumorchid4')

ggplot(alphadiv, aes(x = captivity_status, y = PD, fill=captivity_status))  + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none") +
  ylab("Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = ape_palette)
ggsave('analyses/figures/FigSX_alpha_div_16S_PD.pdf')


#BETA DIVERSITY
#Unweighted unifrac
description_color_scale <- levels(as.factor(metadata$Description))
(description_color_scale <- as.vector(recode(description_color_scale, wild_gorilla = "darkgreen", 
                                         wild_chimp = "darkorange2",
                                         wild_bonobo = "red2", 
                                         western_human = "blue",
                                         non_western_human="skyblue2",
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

library(cowplot)
plot_grid(pcoa_unweighted_unifrac,pcoa_weighted_unifrac, labels = c('A', 'B'), ncol=1)
ggsave('analyses/figures/FigSX_PCoA_unweighted_weighted_unifrac.pdf',width=6)

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

library(PMCMRplus)
library(reshape2)
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

get_kw("Oribacterium")

genera = unique(physeq_gen_abund$Genus)
#calculate for all genera
df = tibble()
for (genus in unique(physeq_gen_abund$Genus)) {
  print(genus)
  res <- get_kw(genus)
  df = df %>% bind_rows(res)
}

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

#adjust pvalues for multiple comparisons
df$p.adj <- round(p.adjust(df$value,method = "fdr",n = length(df$value)),3)
length(unique(df$Genus))

#reshape dataframe so each row is one genus
df <- df %>% select(Genus,comparison,p.adj) %>% spread(comparison,p.adj)

#add mean abundance of each captivity status group
df$wild_mean <- mapply(get_mean, df$Genus, 'wild')
df$captive_mean <- mapply(get_mean, df$Genus, 'captive')
df$nwh_mean <- mapply(get_mean, df$Genus, 'non_western_human')
df$wh_mean <- mapply(get_mean, df$Genus, 'western_human')

#add taxonomy 
df <- physeq_gen_abund %>% 
  select(Phylum,Class,Order,Family,Genus) %>% 
  unique() %>% 
  left_join(df,by='Genus')
dim(df) 
write.table(df,file='analyses/tables/Table_differential_abundance_genera_16S.txt',sep='\t',row.names = F)

#genera that differ significantly in at least one pairwise comparison 
df_sig <- df %>% filter(non_western_humanvs.captive<.05|western_humanvs.captive<.05|western_humanvs.non_western_human<.05|
                          wildvs.captive<.05|wildvs.non_western_human<.05|wildvs.western_human<.05) 

#genera that differ significantly in at least one pairwise comparison  and reach 1% abundance in at least one group
df_sig_over01 <- df %>% filter(non_western_humanvs.captive<.05|western_humanvs.captive<.05|western_humanvs.non_western_human<.05|
                                 wildvs.captive<.05|wildvs.non_western_human<.05|wildvs.western_human<.05) %>%
  filter(wild_mean>.01|captive_mean>.01|nwh_mean>.01|wh_mean>.01)

#genera that reach 5% abundance in at least one group
df_over05 <-  df %>% filter(wild_mean>.05|captive_mean>.05|nwh_mean>.05|wh_mean>.05)
select_gen <- physeq_gen_abund %>% filter(Genus %in% df_over05$Genus) 
select_gen$captivity_status <-ordered(select_gen$captivity_status, levels = c('wild','captive','non_western_human','western_human'))
select_gen$Genus <-ordered(select_gen$Genus, levels = c("Treponema_2","Prevotella_9","Ruminococcaceae_UCG-005",
                                                        "Bacteroides","Blautia","Faecalibacterium",
                                                        "Acinetobacter","Coriobacteriaceae_UCG-003","Flexilinea"))
ape_palette = c("darkorange1",'deepskyblue',"goldenrod1",'mediumorchid4')
df_over05$Genus
select_gen$captivity_status <- as.vector(recode(select_gen$captivity_status, wild = 'wild ape', 
                                                       captive = 'captive ape',
                                                       western_human = "western human",
                                                       non_western_human="non western human"))
new_order = c("wild ape","captive ape","non western human","western human")
select_gen$captivity_status <- factor(select_gen$captivity_status, levels = new_order)

ggplot(select_gen)  + 
  geom_boxplot(aes(x = captivity_status, y = relative_abundance, fill=captivity_status),outlier.shape = NA) + 
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = ape_palette) +
  facet_wrap("Genus",scale="free",nrow=3) +
  theme(axis.title.x = element_blank()) +
  ylab("Relative abundance") 
ggsave('analyses/figures/FigSX_relative_abundance_common_genera_16S.pdf',width=9)
