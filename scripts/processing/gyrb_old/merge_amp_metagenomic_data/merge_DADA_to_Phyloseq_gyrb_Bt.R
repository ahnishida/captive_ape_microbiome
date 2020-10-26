library(dada2); packageVersion("dada2")
library(phyloseq)
library(tidyverse)
library(phytools)
library(seqinr)
library("genefilter")

setwd('/Volumes/AHN/captive_ape_microbiome')
outpath='results/gyrb_all/Bt'
dir.create(outpath, recursive=TRUE)

#import data from separate gyrb amplicon runs and human metagenomic data
print('import gyrb_moeller_wild dataset')
gyrb_moeller_wild_path <- "results/gyrb_moeller_wild/DADA2/Bt" 
gyrb_moeller_wild.seqtab.nochim <- readRDS(file.path(gyrb_moeller_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_moeller_wild.seqtab.nochim )
print('import gyrb_nishida_captive_wild dataset')
gyrb_nishida_captive_wild_path <- "results/gyrb_nishida_captive_wild/DADA2/Bt" 
gyrb_nishida_captive_wild.seqtab.nochim <- readRDS(file.path(gyrb_nishida_captive_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_nishida_captive_wild.seqtab.nochim)
print('import gyrb_nishida_projectchimps dataset')
gyrb_nishida_projectchimps_path <- "results/gyrb_nishida_projectchimps/DADA2/Bt" 
gyrb_nishida_projectchimps.seqtab.nochim <- readRDS(file.path(gyrb_nishida_projectchimps_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_nishida_projectchimps.seqtab.nochim )
print("import human_metagenomic_data")
metagenomic_path <- "results/gyrb_metagenomic/gyrb_metagenomic_seqtab.txt"
metagenomic_human_reads <- as.matrix(read.table(metagenomic_path,header=T,
                                                row.names=1, check.names=F, sep="\t")) #Opens ASV table
rownames(metagenomic_human_reads)
colnames(metagenomic_human_reads)
rowSums(metagenomic_human_reads)
#merge ASV tables
gyrb.asv <- mergeSequenceTables(gyrb_moeller_wild.seqtab.nochim,
                                gyrb_nishida_captive_wild.seqtab.nochim,
                                gyrb_nishida_projectchimps.seqtab.nochim,
                                metagenomic_human_reads)
saveRDS(gyrb.asv,file.path(outpath,"gyrb.asv.rds"))
print("raw asv table dimensions")
print(dim(gyrb.asv))

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(gyrb.asv)
asv_headers <- vector(dim(gyrb.asv)[2], mode="character")
for (i in 1:dim(gyrb.asv)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(outpath,"ASVs_all.fasta"))

#count table:
asv_tab <- t(gyrb.asv)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(outpath,"ASVs_all_counts.tsv"), sep="\t", quote=F, col.names=NA)
dim(asv_tab)

#
tax_table <- read.table(paste0(outfolder,'/ASVs_taxonomy.txt'),header = TRUE)
tax_table_ASVs <- tax_table$ASV

#output select ASV fasta
(select_asv_indices = which(asv_headers %in% paste0(">", tax_table_ASVs)))
(select_asv_headers = asv_headers[select_asv_indices])
(select_asv_seqs = asv_seqs[select_asv_indices])
(select_asv_fasta <-  c(rbind(select_asv_headers, select_asv_seqs)))
write(select_asv_fasta, file.path(outpath,"ASVs_filtered.fasta"))

#assign taxonomy using GTDBTK database
prot_db='results/gyrb_bt_gtdbtk_ref/blast_db/gtdbtk_gyrb_amplicon.faa'
nucl_db='results/gyrb_bt_gtdbtk_ref/blast_db/gtdbtk_gyrb_amplicon.fasta'
asv_fasta_path=file.path(outpath,"ASVs_all.fasta")
asv_faa_path=file.path(outpath,"ASVs_all.faa")
print('assigning taxonomy')
system(paste0('mkdir ',outpath,'/assigned_taxonomy'))
system(paste0('blastn -query ',asv_fasta_path,
             ' -db ',nucl_db,
             ' -outfmt "7 qseqid sacc pident qlen length evalue bitscore salltitles sseq" ',
             ' -max_target_seqs 1 ',
             ' -out ',outpath,"/assigned_taxonomy/ASVs_blastn.txt"))

system(paste0('transeq -frame 2 -sequence ',asv_fasta_path,' -outseq ',asv_faa_path))
system(paste0('blastp -query ',asv_faa_path,
              ' -db ',prot_db,
              ' -outfmt "7 qseqid sacc pident qlen length evalue bitscore salltitles sseq qseq" ',
              ' -max_target_seqs 1 ',
              ' -out ',outpath,"/assigned_taxonomy/ASVs_blastp.txt"))


ASV <- sub(">", "", asv_headers)
ASV <- data.frame(ASV)
pident_cutoff <- 50
length_cutoff <- .8

asv_faa <- read.fasta(asv_faa_path,as.string = TRUE,seqonly=TRUE)
asv_names <- sub(">", "", asv_headers)
(ASVs_w_stop_codons <- asv_names[str_detect(asv_faa,'\\*')])

blastp <- read.table(file.path(outpath,"assigned_taxonomy/ASVs_blastp.txt")) 
colnames(blastp) <- c('ASV','sacc','pident','qlen','length','evalue','bitscore','genome','salltitles','species_id','sseq','qseq')
blastp$ASV <- str_sub(blastp$ASV, 1, str_length(blastp$ASV)-2)
blastp <- blastp %>% separate(salltitles,sep=';',into=c('Domain','Phylum','Class','Order',
                                                        'Family','Genus','Species')) %>%
  filter(pident>pident_cutoff & length>(83*length_cutoff))  %>%
  filter(!ASV %in% ASVs_w_stop_codons) %>%
  filter(Order=='o__Bacteroidales') %>%
  select('ASV','Phylum','Class','Order','Family')

blastn <- read.table(file.path(outpath,"assigned_taxonomy/ASVs_blastn.txt"))
colnames(blastn) <- c('ASV','sacc','pident','qlen','length','evalue','bitscore','genome','salltitles','species_id','sseq')
blastn <- blastn %>% separate(salltitles,sep=';',into=c('Domain','Phylum','Class','Order',
                                                        'Family','Genus','Species')) %>% 
                      filter(pident>pident_cutoff & length>(250*length_cutoff)) %>% 
                      mutate(Genus_pident_cutoff = ifelse(pident>80,Genus,'Unassigned')) %>%
                      mutate(Species_pident_cutoff = ifelse(pident>95,Species,'Unassigned')) %>%
                      select('ASV','Genus','Species')

tax_table <- ASV %>% 
  left_join(blastp, by='ASV') %>% 
  left_join(blastn, by='ASV') %>%
  mutate_each(funs(replace(., which(is.na(.)), 'Unassigned')))

tax_table_filt <- tax_table %>% filter(Family != 'Unassigned')
write.table(tax_table_filt, file=file.path(outpath,"assigned_taxonomy/ASVs_filtered_taxonomy.txt"),sep='\t',quote=F,row.names=F)
taxa = read.csv(file.path(outpath,"assigned_taxonomy/ASVs_filtered_taxonomy.txt"),sep='\t')
taxa = taxa %>% column_to_rownames(var='ASV')

#output select ASV fasta
(select_gen_ASVs = rownames(taxa)) #get list of ASVs belonging to select taxa 
(select_asv_indices = which(asv_headers %in% paste0(">", select_gen_ASVs)))
(select_asv_headers = asv_headers[select_asv_indices])
(select_asv_seqs = asv_seqs[select_asv_indices])
(select_asv_fasta <-  c(rbind(select_asv_headers, select_asv_seqs)))
write(select_asv_fasta, file.path(outpath,"ASVs_filtered.fasta"))

#output select ASV count table
select_asv_tab <- asv_tab[row.names(asv_tab) %in% select_gen_ASVs, ]
print('Bacteroidales ASV table')
dim(select_asv_tab)
write.table(select_asv_tab, file.path(outpath,"ASVs_filtered_counts.tsv"),
            sep="\t", quote=F, col.names=NA)

#generate phylogeny
system("mkdir results/gyrb_all/Bt/phylogeny")
system("sed 's/[; ]//g' results/gyrb_bt_gtdbtk_ref/alignment/gtdbtk_gyrb_Bt_amplicon.faa > results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt_amplicon.faa")
system("sed 's/[; ]//g' results/gyrb_bt_gtdbtk_ref/alignment/gtdbtk_gyrb_Bt_amplicon.fasta > results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt_amplicon.fasta")
system("sed 's/[; ]//g' results/gyrb_bt_gtdbtk_ref/alignment/gtdbtk_gyrb_Bt.faa > results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt.faa")
system("sed 's/[; ]//g' results/gyrb_bt_gtdbtk_ref/alignment/gtdbtk_gyrb_Bt.fasta > results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt.fasta")


system("transeq -frame 2 -sequence results/gyrb_all/Bt/ASVs_filtered.fasta -outseq results/gyrb_all/Bt/ASVs_filtered.faa")
system("cat results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt_amplicon.faa results/gyrb_all/Bt/ASVs_filtered.faa >  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.faa ")
system("cat results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt_amplicon.fasta results/gyrb_all/Bt/ASVs_filtered.fasta >  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.fasta ")
system("cat results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt.faa results/gyrb_all/Bt/ASVs_filtered.faa >  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.faa ")
system("cat results/gyrb_all/Bt/phylogeny/gtdbtk_gyrb_Bt.fasta results/gyrb_all/Bt/ASVs_filtered.fasta >  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.fasta ")

system("mafft results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.faa > results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.faa.aln")
system("tranalign -asequence  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.fasta -bsequence results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.faa.aln -outseq results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.fasta.aln")
system("mafft results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.faa > results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.faa.aln")
system("tranalign -asequence  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.fasta -bsequence results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.faa.aln -outseq results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.fasta.aln")

system("./FastTree-2.1.9 -nt -gtr <  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.fasta.aln >  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref.fasta.aln.tree")
system("./FastTree-2.1.9 -nt -gtr <  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.fasta.aln >  results/gyrb_all/Bt/phylogeny/ASVs_filtered_ref_fulllength.fasta.aln.tree")


#import into phyloseq
TAX = tax_table(as.matrix(taxa))
OTU = otu_table(select_asv_tab,taxa_are_rows = TRUE)
TREE =  read.tree(file.path(outpath,"phylogeny/ASVs_filtered_ref.fasta.aln.tree"))


#format sample metadata
all_samples_metadata_file = "metadata/metadata_all_samples_gyrb.txt"
all_samples_metadata <- read.csv(all_samples_metadata_file,sep='\t')
all_samples_metadata$X.SampleID
metadata = all_samples_metadata %>%  #filter metadata based on whether samples are in OTU table
  filter(X.SampleID %in% sample_names(OTU)) %>%
  column_to_rownames('X.SampleID')
write.table(metadata,file = "metadata/metadata_Bt_samples_gyrb.txt",sep='\t',quote=FALSE)
setdiff(sample_names(OTU),all_samples_metadata$X.SampleID) #verify no samples missing from metadata
SAM = sample_data(metadata)
metadata_summary <- metadata %>% group_by(Description,site) %>% tally()

#root phylogeny
ps <- phyloseq(OTU,SAM,TAX,TREE)
plot_tree(ps, color="Family")
outgroup<-subset_taxa(ps,Family=='f__UBA932')
outgroup_taxa <- findMRCA(phy_tree(ps),taxa_names(outgroup))
phy_tree(ps) <- root(phy_tree(ps),node=outgroup_taxa)
plot_tree(ps, color="Family")
plot_tree(ps, color="Genus")
plot_tree(ps, color="Description")
saveRDS(ps,file.path(outpath,'Bacteroidales_physeq.rds'))

#view just Bacteroidaceae
Bacteroidaceae <- subset_taxa(ps,Genus=='g__Prevotella'|Genus=='g__Bacteroides')
(Bacteroidaceae_node <- findMRCA(phy_tree(ps),taxa_names(Bacteroidaceae)))
Bacteroidaceae_tree <- extract.clade(phy_tree(ps),Bacteroidaceae_node)
Bacteroidaceae <- phyloseq(OTU,SAM,TAX,Bacteroidaceae_tree)

(Bacteroidaceae <- subset_taxa(ps,Family=='f__Bacteroidaceae'))
unique(as.data.frame(tax_table(Bacteroidaceae))$Genus)
plot_tree(Bacteroidaceae, color="Genus")
saveRDS(ps,file.path(outpath,'Bacteroidaceae_physeq.rds'))

# Histogram of sample read counts
sample_sum_df <- data.frame(sum = sample_sums(Bacteroidaceae)) %>% rownames_to_column('SampleID')
metadata <- metadata %>% rownames_to_column('SampleID')
metadata$SampleID <- as.character(metadata$SampleID)
sample_sum_df  = sample_sum_df %>% full_join(metadata,by='SampleID')
sample_sum_df$dataset <- as.character(sample_sum_df$dataset)

ggplot(sample_sum_df, aes(x = sum,color=dataset,fill=dataset)) + 
  geom_histogram(color = "black",binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#relative abundance filtering
FSr  = transform_sample_counts(Bacteroidaceae, function(x) x / sum(x) )
flist<-filterfun(kOverA(2,.0))
(k2overa0 =filter_taxa(FSr,flist,TRUE))
(Bacteroidaceae_k2overa0<- prune_taxa(taxa_names(k2overa0),Bacteroidaceae))
plot_tree(Bacteroidaceae_k2overa0,color='Genus')
saveRDS(Bacteroidaceae_k2overa0,file.path(outpath,'Bacteroidaceae_physeq_k2overa0.rds'))

#Prevotella super clade, 
#This includes  wild ape ASVs that were not assigned at genus level bc of no close sequence
Prevotella_taxa <- subset_taxa(Bacteroidaceae_k2overa0,Genus=='g__Prevotella'|Genus=='g__UBA4334')
(Prevotella_node <- findMRCA(phy_tree(Bacteroidaceae_k2overa0),taxa_names(Prevotella_taxa)))
Prevotella_tree <- extract.clade(phy_tree(Bacteroidaceae_k2overa0),Prevotella_node)
Prevotella <- phyloseq(OTU,SAM,TAX,Prevotella_tree)
unique(as.data.frame(tax_table(Prevotella))$Genus)
plot_tree(Prevotella, color="Genus")
saveRDS(Prevotella ,file.path(outpath,'Prevotella_phyloseq.rds'))


pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

#Generating a paraphyletic Bacteroides tree
(Bacteroides <- subset_taxa(Bacteroidaceae_k2overa0,Genus=='g__Bacteroides'))
(Bacteroides_node <- findMRCA(phy_tree(Bacteroidaceae_k2overa0),taxa_names(Bacteroidaceae_k2overa0)))
Bacteroides_tree <- extract.clade(phy_tree(Bacteroidaceae_k2overa0),Bacteroides_node)
(Bacteroides <- phyloseq(OTU,SAM,TAX,Bacteroides_tree))
Bacteroides_collapse <- merge_samples(Bacteroides,group = "Description")
plot_tree(Bacteroides_collapse, color="Genus")

#filter out the clade that contains Prevotella and other genera
(Prevotella_taxa <- subset_taxa(Bacteroides,Genus=='g__Prevotella'|Genus=='g__UBA4334'))
(Other_genera_node <- findMRCA(phy_tree(Bacteroides),taxa_names(Other_genera)))
Other_genera_names <- extract.clade(phy_tree(Bacteroides),Other_genera_node)
(Bacteroides_para <- pop_taxa(Bacteroides,Other_genera_names$tip.label))
plot_tree(Bacteroides_para, color="Genus")
saveRDS(Bacteroides_para, file.path(outpath,'Bacteroides_phyloseq.rds'))


