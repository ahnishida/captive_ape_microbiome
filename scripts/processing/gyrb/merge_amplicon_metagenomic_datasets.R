library(dada2); packageVersion("dada2")
library(phyloseq)
library(DECIPHER)
library(tidyverse)
library(phytools)
library(seqinr)
library(stringr)
#library("genefilter")
setwd('/Volumes/AHN/captive_ape_microbiome')

# this merges amplicon datasets and filters out reads that aren;t assign. 
# filtered ASVs are then merged with the gtdbtk reference fasta and used to blastn query metagenomic datasets

#load assign taxonomy script
source('scripts/processing/gyrb/idTaxa.R')

#define inputs/outputs
gyrb_ref_faa = 'ref_seqs/gtdbtk_gyrb.faa'
outpath='results/gyrb/processing/gyrb_amp_meta_datasets'
dir.create(outpath, recursive=TRUE)

#import data from separate gyrb amplicon runs and human metagenomic data
print('import gyrb_moeller_wild dataset')
gyrb_moeller_wild_path <- "results/gyrb/processing/gyrb_moeller_wild/DADA2" 
gyrb_moeller_wild.seqtab.nochim <- readRDS(file.path(gyrb_moeller_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_moeller_wild.seqtab.nochim )
print('import gyrb_nishida_captive_wild dataset')
gyrb_nishida_captive_wild_path <- "results/gyrb/processing/gyrb_nishida_captive_wild/DADA2" 
gyrb_nishida_captive_wild.seqtab.nochim <- readRDS(file.path(gyrb_nishida_captive_wild_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_nishida_captive_wild.seqtab.nochim)
print('import gyrb_nishida_projectchimps dataset')
gyrb_nishida_projectchimps_path <- "results/gyrb/processing/gyrb_nishida_projectchimps/DADA2" 
gyrb_nishida_projectchimps.seqtab.nochim <- readRDS(file.path(gyrb_nishida_projectchimps_path,"seqtab.nochim.rds")) #Opens ASV table
rownames(gyrb_nishida_projectchimps.seqtab.nochim )
print("import human_metagenomic_data")
metagenomic_path <- "results/gyrb/processing/metagenomic_samples/final_outputs/gyrb_metagenomic_seqtab.txt"
metagenomic_human_reads <- as.matrix(read.table(metagenomic_path,header=T,
                                                row.names=1, check.names=F, sep="\t")) #Opens ASV table
dim(metagenomic_human_reads)
saveRDS(seqtab,file.path(outpath,'seqtab.rds'))
seqtab = readRDS(file.path(outpath,'seqtab.rds'))
#merge ASV tables
seqtab <- mergeSequenceTables(gyrb_moeller_wild.seqtab.nochim,
                                gyrb_nishida_captive_wild.seqtab.nochim,
                                gyrb_nishida_projectchimps.seqtab.nochim,
                                metagenomic_human_reads)
saveRDS(seqtab,file.path(outpath,'seqtab.rds'))
seqtab = readRDS(file.path(outpath,'seqtab.rds'))

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(outpath,"ASVs_all.fasta"))
rownames(seqtab)

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)

# metadata
all_samples_metadata_file = "metadata/metadata_gyrb_amp_meta.txt"
metadata <- read.csv(all_samples_metadata_file,sep='\t')
metadata = metadata %>%  #filter metadata based on whether samples are in OTU table
  filter(X.SampleID %in% colnames(asv_tab)) %>%
  column_to_rownames('X.SampleID')
setdiff(colnames(asv_tab),row.names(metadata)) #samples in otu table that aren't in metadata

#filter ASVs
asv_fasta <- readDNAStringSet(file.path(outpath,"ASVs_all.fasta"),format = 'fasta')
asv_faa <- Biostrings::translate(DNAStringSet(asv_fasta,start =2))
Biostrings::writeXStringSet(asv_faa,file=file.path(outpath,"ASVs_all.faa")) 
system(paste('./scripts/processing/gyrb/blastp_filter_ASVs.sh',
      file.path(outpath,"ASVs_all.faa"),
      gyrb_ref_faa, 
      file.path(outpath,"ASVs_all_blastp.txt"),sep=' '))
ASV_blastp_res <- read.table(file.path(outpath,"ASVs_all_blastp.txt"))
colnames(ASV_blastp_res) <- c('ASV','sacc','pident','qlen','length','evalue','bitscore','genome',
                              'salltitles','species','sseq','qseq','sstart','send')

#filter ASV based on alignment length and percent identity
ASVs_w_stop_codons <- names(asv_faa[str_detect(asv_faa,'\\*')]) #identity asvs that when translated have stop codons
ASVs_filtered_blast <- ASV_blastp_res %>%
  filter(length>(83*.90)&length<(83*1.10))  %>%
  filter(pident > 80)  %>%
  filter(!ASV %in% ASVs_w_stop_codons)
#there are a few outliers not aligning to the same region as the majority of ASVs
hist(ASVs_filtered_blast$sstart)
hist(ASVs_filtered_blast$send)
#filter ASVs based on where they hit to the alignment 
ASVs_filtered_blast <-ASVs_filtered_blast %>% 
  filter(sstart<150) %>% 
  filter(send<220)
hist(ASVs_filtered_blast$sstart)
hist(ASVs_filtered_blast$send)

#output ASVs filtered based on blastp results 
ASVs_filtered_blast_fasta <- asv_fasta[names(asv_fasta) %in% ASVs_filtered_blast$ASV]
names(ASVs_filtered_blast_fasta)
Biostrings::writeXStringSet(ASVs_filtered_blast_fasta ,file=file.path(outpath,"ASVs_filtered_blast.fasta")) 

#assign taxonomy
TAX_table <- assign_taxonomy_w_idTAXA(file.path(outpath,"ASVs_all.faa"),gyrb_ref_faa)
data.frame(TAX_table) %>% group_by(Family) %>% tally()
outgroup <- data.frame(TAX_table) %>% filter(Family == 'f__F082') 
F082_ASVs <- rownames(outgroup) 

# import into phyloseq
(ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
                sample_data(metadata), 
                tax_table(as.matrix(TAX_table))))
saveRDS(ps,file.path(outpath,'ps.rds'))
ps = readRDS(file.path(outpath,'ps.rds'))

#write physeq
ps <- prune_taxa(ASVs_filtered_blast$ASV,ps)
ps_Bacteroidales <- subset_taxa(ps, Order == 'o__Bacteroidales')
ps_Bacteroidales <- prune_samples(sample_sums(ps_Bacteroidales)>0, ps_Bacteroidales)
saveRDS(ps_Bacteroidales,file.path(outpath,"ps_Bacteroidales.rds"))
ps_Bacteroidales <- readRDS(file.path(outpath,"ps_Bacteroidales.rds"))

#write fasta
asv_fasta_filt <- asv_fasta[names(asv_fasta) %in% taxa_names(ps_Bacteroidales)]
Biostrings::writeXStringSet(asv_fasta_filt,file=file.path(outpath,"ASVs_filtered.fasta")) 
asv_faa_filt <- asv_faa[names(asv_faa) %in% taxa_names(ps_Bacteroidales)]
Biostrings::writeXStringSet(asv_faa_filt,file=file.path(outpath,"ASVs_filtered.faa")) 

#select ref seqs to include with ASVs
gyrb_ref_fasta = 'ref_seqs/gtdbtk_gyrb.fasta'
gyrb_ref <- readDNAStringSet(gyrb_ref_fasta,format = 'fasta')
bacteroidales_ref <- gyrb_ref[str_detect(names(gyrb_ref),'o__Bacteroidales')]
Biostrings::writeXStringSet(bacteroidales_ref ,file=file.path(outpath,"bacteroidales_ref.fasta")) 

#phylogeny 
folder <- 'phylogeny'
dir.create(file.path(outpath,folder), recursive=TRUE)

trimAlign <- function(aln_in,aln_out){
  #trim alignment to start and end of ASVs
  aln <- Biostrings::readDNAStringSet(aln_in)
  ASV_1 <- DNAString(as.character(aln['ASV_1']))
  ASV_seq_no_gaps <- as.character(RemoveGaps(aln['ASV_1'],removeGaps = "all"))
  first10nt <- stringr::str_sub(ASV_seq_no_gaps,1,10)
  last10nt <- stringr::str_sub(ASV_seq_no_gaps,-10,-1)
  s <- start(matchPattern(first10nt, ASV_1))
  e <-end(matchPattern(last10nt, ASV_1))
  alnTrim <- DNAStringSet(aln,start=s,end=e)
  seqlengths = width(RemoveGaps(alnTrim,
                                removeGaps = "all",
                                processors = 1))
  alnTrimFilt <- alnTrim[seqlengths > 250*.95]
  alnTrimFilt <- RemoveGaps(alnTrimFilt,
                            removeGaps = "common")
  Biostrings::writeXStringSet(alnTrimFilt,file=aln_out)
}

#no ref taxa
asv_fasta_filt <- Biostrings::readDNAStringSet(file.path(outpath,"ASVs_filtered.fasta"))
asv_aln <- AlignTranslation(DNAStringSet(asv_fasta_filt,start=2))
Biostrings::writeXStringSet(asv_aln,file.path(outpath,folder,'ASVs_filtered_aln.fasta'))
trimAlign(file.path(outpath,folder,'ASVs_filtered_aln.fasta'),
          file.path(outpath,folder,"ASVs_filtered_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_aln_trim.tre")))
asv_tree <- ape::read.tree(file.path(outpath,folder,"ASVs_filtered_aln_trim.tre"))
mrca <- findMRCA(asv_tree,F082_ASVs)
asv_tree_rooted <- reroot(asv_tree,mrca)
ape::write.tree(asv_tree_rooted,"tmp.tre")
system(paste0("sed 's/Root;/;/g'  tmp.tre > ",
              file.path(outpath,folder,"ASVs_filtered_aln_trim_rooted.tre")))
ps_Bacteroidales_asvTree <- merge_phyloseq(ps_Bacteroidales,asv_tree_rooted)
saveRDS(ps_Bacteroidales_asvTree,file.path(outpath,"ps_Bacteroidales_asvTree.rds"))

#ref taxa
asv_ref <- c(bacteroidales_ref,DNAStringSet(asv_fasta_filt,start=2))
Biostrings::writeXStringSet(asv_ref,file=file.path(outpath,folder,"ASVs_filtered_ref.fasta")) 
asv_ref_aln <- AlignTranslation(asv_ref)  
Biostrings::writeXStringSet(asv_ref_aln,file=file.path(outpath,folder,"ASVs_filtered_ref_aln.fasta")) 
trimAlign(file.path(outpath,folder,'ASVs_filtered_ref_aln.fasta'),
          file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.fasta"),
              ' > ',
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre")))
asv_ref_tree <- ape::read.tree(file.path(outpath,folder,"ASVs_filtered_ref_aln_trim.tre"))
asv_tips <- asv_ref_tree$tip.label[str_detect(asv_ref_tree$tip.label,'ASV')]
asv_ref_tree <- keep.tip(asv_ref_tree,asv_tips)
mrca <- findMRCA(asv_ref_tree,F082_ASVs)
asv_ref_tree_rooted <- reroot(asv_ref_tree,mrca)
ape::write.tree(asv_ref_tree_rooted,"tmp.tre")
system(paste0("sed 's/Root;/;/g'  tmp.tre > ",
              file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre")))
ps_Bacteroidales_asvRefTree <- merge_phyloseq(ps_Bacteroidales,asv_ref_tree_rooted)
saveRDS(ps_Bacteroidales_asvRefTree,file.path(outpath,"ps_Bacteroidales_asvRefTree.rds"))

#any diff b/n physeqs 
setdiff(taxa_names(ps_Bacteroidales_asvTree),taxa_names(ps_Bacteroidales_asvRefTree))
setdiff(taxa_names(ps_Bacteroidales_asvRefTree),taxa_names(ps_Bacteroidales_asvTree))

#write to inputs folder
write_to_output_folder <- function(physeq,outdir,all_ASVs_fasta){
  dir.create(outdir, recursive=TRUE)
  #write tax table
  tax_table_tab <- as.data.frame(tax_table(physeq)) %>% 
    rownames_to_column(var='ASV')
    write.table(tax_table_tab,file.path(outdir,'physeq_Bacteroidales_taxonomy.txt'),quote=F,row.names=F,sep='\t')
  metadata_tab <- as.data.frame(sample_data(physeq)) 
    metadata_tab$X.SampleID <- rownames(metadata_tab)
    write.table(metadata_tab,file.path(outdir,'physeq_metadata_passing_samples.txt'),quote=F,row.names=F,sep='\t')
  otu_table <- as.data.frame(otu_table(physeq)) 
    write.table(otu_table,file.path(outdir,'physeq_Bacteroidales_asv_tab.txt'),quote=F,row.names=T,col.names=NA,sep='\t')
  fasta <- Biostrings::readDNAStringSet(all_ASVs_fasta)
    fasta_physeq <- fasta[names(fasta) %in% taxa_names(physeq)]
    Biostrings::writeXStringSet(fasta_physeq,file=file.path(outdir,"physeq_Bacteroidales_asv.fasta")) 
}

write_to_output_folder(ps_Bacteroidales_asvTree,
                       'results/gyrb/inputs/ps_Bacteroidales_asvTree',
                       file.path(outpath,"ASVs_all.fasta"))
write_to_output_folder(ps_Bacteroidales_asvRefTree,
                       'results/gyrb/inputs/ps_Bacteroidales_asvRefTree',
                       file.path(outpath,"ASVs_all.fasta"))
#copy over phylogenies
system(paste0('cp ',file.path(outpath,folder,"ASVs_filtered_ref_aln_trim_rooted.tre"),
' results/gyrb/inputs/ps_Bacteroidales_asvRefTree/physeq_Bacteroidales.tree'))
system(paste0('cp ',file.path(outpath,folder,"ASVs_filtered_aln_trim_rooted.tre"),
              ' results/gyrb/inputs/ps_Bacteroidales_asvTree/physeq_Bacteroidales.tree'))

#write moeller codiv seq files
system('cp results/gyrb/processing/moeller_sup/moeller_codiv* results/gyrb/inputs/ps_Bacteroidales_asvTree/')
system('cp results/gyrb/processing/moeller_sup/moeller_codiv* results/gyrb/inputs/ps_Bacteroidales_asvRefTree/')
