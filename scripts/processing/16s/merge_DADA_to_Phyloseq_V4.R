library(dada2); packageVersion("dada2")
library("genefilter")
library(tidyverse)
library(phyloseq)
library(seqinr)
library(ape)
library(phytools)

setwd('/Volumes/AHN/captive_ape_microbiome/results/16s')

dir.create("processing/16s__merged", recursive = TRUE)
dir.create("inputs")

# read in 16s datasets ASV tables, trim reads to 250 character so they will all be same length as gor bon dataset
raymann_ASV <- readRDS("processing/16s_raymann_captive/DADA2/seqtab.nochim.rds") #Opens ASV table
hist(nchar(colnames(raymann_ASV  )))
colnames(raymann_ASV) <- substr(colnames(raymann_ASV),start=1, stop=250)
rowSums(raymann_ASV)

clayton_ASV <- readRDS("processing/16s_clayton_captive/DADA2/seqtab.nochim.rds") #Opens ASV table
hist(nchar(colnames(clayton_ASV)))
colnames(clayton_ASV) <- substr(colnames(clayton_ASV),start=1, stop=250)
rowSums(clayton_ASV)

moeller_chimp_ASV <- readRDS("processing/16s_moeller_wild/chimp/DADA2/seqtab.nochim.rds") #Opens ASV table
moeller_chimp_ASV <- moeller_chimp_ASV[,colnames(moeller_chimp_ASV)[nchar(colnames(moeller_chimp_ASV))>=250]]
hist(nchar(colnames(moeller_chimp_ASV)))
colnames(moeller_chimp_ASV) <- substr(colnames(moeller_chimp_ASV),start=1, stop=250)
rowSums(moeller_chimp_ASV)

moeller_gor_bon_ASV <- readRDS("processing/16s_moeller_wild/gor_bon/DADA2/seqtab.nochim.rds") #Opens ASV table
rowSums(moeller_gor_bon_ASV)

nishida_captive_ASV <- readRDS("processing/16s_nishida_captive/DADA2/seqtab.nochim.rds") #Opens ASV table
colnames(nishida_captive_ASV) <- substr(colnames(nishida_captive_ASV),start=1, stop=250)
rowSums(nishida_captive_ASV)

nishida_projectchimps_ASV <- readRDS("processing/16s_nishida_projectchimps/DADA2/seqtab.nochim.rds") #Opens ASV table
colnames(nishida_projectchimps_ASV) <- substr(colnames(nishida_projectchimps_ASV),start=1, stop=250)
rowSums(nishida_projectchimps_ASV)

vangay_ASV <- readRDS("processing/16s_vangay_2018/DADA2/seqtab.nochim.rds") #Opens ASV table
colnames(vangay_ASV) <- substr(colnames(vangay_ASV),start=1, stop=250)
dim(vangay_ASV)

smits_ASV <- readRDS("processing/16s_smits_2017/DADA2/seqtab.nochim.rds") #Opens ASV table
colnames(smits_ASV) <- substr(colnames(smits_ASV),start=1, stop=250)
rowSums(smits_ASV)

goodrich_ASV <- readRDS("processing/16s_goodrich/DADA2/seqtab.nochim.rds") #Opens ASV table
colnames(goodrich_ASV) <- substr(colnames(goodrich_ASV),start=1, stop=250)
rowSums(goodrich_ASV)

# merge ASV table 
seqtab <- mergeSequenceTables(raymann_ASV,
                              clayton_ASV,
                              moeller_chimp_ASV,
                              moeller_gor_bon_ASV,
                              nishida_captive_ASV,
                              nishida_projectchimps_ASV,
                              vangay_ASV,
                              smits_ASV,
                              goodrich_ASV)
dim(seqtab) #

# assign taxonomy
TAX <- assignTaxonomy(seqtab, "../../ref_seqs/silva_nr_v132_train_set.fa", multithread=TRUE)
saveRDS(TAX,"processing/16s__merged/tax_table.rds")
TAX = readRDS("processing/16s__merged/tax_table.rds")
dim(TAX)
TAX.print <- TAX # Removing sequence rownames for display only
rownames(TAX.print) <- NULL
head(TAX.print)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# correcting sample names so they don't have funny extensions
samples.out <- rownames(seqtab)
samples.out <- sapply(strsplit(basename(samples.out), "_"), `[`, 1)
rownames(seqtab) <- samples.out

# making and writing out a fasta of ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path("processing/16s__merged/ASVs_prefilter.fasta"))

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)

# tax table:
asv_tax <- TAX
row.names(asv_tax) <- sub(">", "", asv_headers)

# metadata
all_samples_metadata_file = "../../metadata/metadata_all_samples_16s.txt"
metadata <- read.csv(all_samples_metadata_file,sep='\t')
metadata = metadata %>%  #filter metadata based on whether samples are in OTU table
  filter(X.SampleID %in% colnames(asv_tab)) %>%
  column_to_rownames('X.SampleID')
setdiff(colnames(asv_tab),row.names(metadata)) #samples in otu table that aren't in metadata
setdiff(row.names(metadata),colnames(asv_tab)) #samples in metadata, not in otu table

# import into phyloseq
(ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
               sample_data(metadata), 
               tax_table(asv_tax)))

# taxon filtering
(ps_filtered <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" |
      Kingdom == "Archaea" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast"))

# relative abundance filtering must be present in 2 individuals at least .5% relative abundance in each
FSr  = transform_sample_counts(ps_filtered, function(x) x / sum(x) )
flist<-filterfun(kOverA(2,.005))
(k2overa.005 =filter_taxa(FSr,flist,TRUE))
(ps_filtered_relabund <- prune_taxa(taxa_names(k2overa.005),ps_filtered))

# rarefaction 10000 
set.seed(1389)
physeq_rare10000 = prune_samples(sample_sums(ps_filtered_relabund)>=10000, ps_filtered_relabund)
(physeq_rare10000 = rarefy_even_depth(physeq_rare10000,trimOTUs = FALSE))

# subsample humans datasets, so that groups have largely equal numbers
metadata <- data.frame(sample_data(physeq_rare10000)) %>% 
  tibble::rownames_to_column(var = 'X.SampleID')
metadata %>% group_by(Description,dataset) %>% tally()
smits_subset = metadata %>%
  filter(dataset == 'smits_2017') %>%  
  sample_n(100)
vangay_subset_wh = metadata %>%
  filter(dataset == 'vangay_2018' & Description == 'industrialized_human') %>%  
  sample_n(70)
vangay_subset_nwh = metadata %>% 
  filter(dataset == 'vangay_2018' & Description == 'non_industrialized_human') 
goodrich_subset = metadata %>%
  filter(dataset == 'goodrich_2016') %>%  
  sample_n(70)
human_metadata_subset = bind_rows(smits_subset,vangay_subset_wh,vangay_subset_nwh,goodrich_subset)
human_metadata_subset %>% group_by(Description,dataset) %>% tally()
other_groups <-  metadata %>%
  filter(Description != 'non_industrialized_human' & Description != 'industrialized_human')  
metadata_subset <- bind_rows(human_metadata_subset,other_groups)
write.table(metadata_subset,sep='\t',quote=F,row.names=F,
            file=file.path('processing/16s__merged/metadata_surviving_samples.txt'))

physeq_rare10000 <- prune_samples(as.vector(metadata_subset$X.SampleID), physeq_rare10000)
taxa_over_zero <- taxa_names(physeq_rare10000)[taxa_sums(physeq_rare10000)>0]
physeq_rare10000 <- prune_taxa(taxa_over_zero,physeq_rare10000)

asv_fasta <- read.fasta(file.path("processing/16s__merged/ASVs_prefilter.fasta"))
asv_filtered_fasta <- asv_fasta[c(which(names(asv_fasta) %in% rownames(otu_table(physeq_rare10000))))]
write.fasta(sequences = asv_filtered_fasta, names = names(asv_filtered_fasta), 
            file.out = "processing/16s__merged/ASVs.fasta")
length(asv_filtered_fasta)
# align seqs and generate phylogeny
system("mafft --auto processing/16s__merged/ASVs.fasta > processing/16s__merged/ASVs_aligned.fasta")
#have to have an executable fasttree to run in Rstudio
system("./../../FastTree-2.1.9 -nt -gtr < processing/16s__merged/ASVs_aligned.fasta > processing/16s__merged/ASVs_aligned.tre")
TREE <- read.tree("processing/16s__merged/ASVs_aligned.tre")
TREE <- midpoint.root(TREE)
TREE$tip.label
#add tree to phyloseq
physeq_rare10000 <- merge_phyloseq(physeq_rare10000,phy_tree(TREE))
saveRDS(physeq_rare10000,"processing/16s__merged/phyloseq_rare10000.rds")

#write to inputs folder
saveRDS(physeq_rare10000,"inputs/phyloseq_rare10000.rds")
tax_table_tab <- as.data.frame(tax_table(physeq_rare10000)) %>% 
  rownames_to_column(var='ASV')
write.table(tax_table_tab,'inputs/ASVs_taxonomy.txt',quote=F,row.names=F,sep='\t')

metadata_tab <- as.data.frame(sample_data(physeq_rare10000)) 
metadata_tab$X.SampleID <- rownames(metadata_tab)
write.table(metadata_tab,'inputs/16S_metadata.txt',quote=F,row.names=F,sep='\t')

otu_table <- as.data.frame(otu_table(physeq_rare10000)) 
write.table(otu_table,'inputs/ASV_tab.txt',quote=F,row.names=T,col.names=NA,sep='\t')

write.tree(phy_tree(physeq_rare10000),'inputs/ASVs.tree')  

