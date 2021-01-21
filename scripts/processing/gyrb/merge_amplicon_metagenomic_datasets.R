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
processing_folder='results/gyrb/processing/gyrb_amp_meta_datasets'
dir.create(processing_folder, recursive=TRUE)

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
#merge ASV tables
seqtab <- mergeSequenceTables(gyrb_moeller_wild.seqtab.nochim,
                                gyrb_nishida_captive_wild.seqtab.nochim,
                                gyrb_nishida_projectchimps.seqtab.nochim,
                                metagenomic_human_reads)
saveRDS(seqtab,file.path(processing_folder,'seqtab.rds'))
seqtab = readRDS(file.path(processing_folder,'seqtab.rds'))

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(processing_folder,"ASVs_all.fasta"))
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

#translate ASVs to prot seq and blast 
asv_fasta <- readDNAStringSet(file.path(processing_folder,"ASVs_all.fasta"),format = 'fasta')
asv_faa <- Biostrings::translate(DNAStringSet(asv_fasta,start =2)) #translate ASVs
Biostrings::writeXStringSet(asv_faa,file=file.path(processing_folder,"ASVs_all.faa")) 
system(paste('./scripts/processing/gyrb/blastp_filter_ASVs.sh', #run blastp of ASVs
      file.path(processing_folder,"ASVs_all.faa"),
      gyrb_ref_faa, 
      file.path(processing_folder,"ASVs_all_blastp.txt"),sep=' '))
ASV_blastp_res <- read.table(file.path(processing_folder,"ASVs_all_blastp.txt"))
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

#assign taxonomy
TAX_table <- assign_taxonomy_w_idTAXA(file.path(processing_folder,"ASVs_all.faa"),gyrb_ref_faa)
data.frame(TAX_table) %>% group_by(Family) %>% tally()
outgroup <- data.frame(TAX_table) %>% filter(Family == 'f__F082') 
F082_ASVs <- rownames(outgroup) 

# import into phyloseq
(ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
                sample_data(metadata), 
                tax_table(as.matrix(TAX_table))))
saveRDS(ps,file.path(processing_folder,'ps.rds'))
ps = readRDS(file.path(processing_folder,'ps.rds'))
ps_filtered <- prune_taxa(ASVs_filtered_blast$ASV,ps) #filter ASVs based on blastp results
ps_Bacteroidales <- subset_taxa(ps_filtered, Order == 'o__Bacteroidales') #filter ASVs based on taxonomy
samples_with_Bacteroidales <- sample_names(ps_Bacteroidales)[sample_sums(ps_Bacteroidales)>0]
ps_Bacteroidales <- prune_samples(samples_with_Bacteroidales, ps_Bacteroidales)

#write fasta
asv_fasta_filt <- asv_fasta[names(asv_fasta) %in% taxa_names(ps_Bacteroidales)]
Biostrings::writeXStringSet(asv_fasta_filt,file=file.path(processing_folder,"ASVs_filtered.fasta")) 
asv_faa_filt <- asv_faa[names(asv_faa) %in% taxa_names(ps_Bacteroidales)]
Biostrings::writeXStringSet(asv_faa_filt,file=file.path(processing_folder,"ASVs_filtered.faa")) 

#select ref seqs to include with ASVs
gyrb_ref_fasta = 'ref_seqs/gtdbtk_gyrb.fasta'
gyrb_ref <- readDNAStringSet(gyrb_ref_fasta,format = 'fasta')
bacteroidales_ref <- gyrb_ref[str_detect(names(gyrb_ref),'o__Bacteroidales')]
Biostrings::writeXStringSet(bacteroidales_ref ,file=file.path(processing_folder,"bacteroidales_ref.fasta")) 

#phylogeny 
dir.create(file.path(processing_folder,'phylogeny'), recursive=TRUE)

trimAlign <- function(aln_in,aln_out){
  #trim alignment to start and end of ASVs
  aln <- Biostrings::readDNAStringSet(aln_in) 
  ASV1 <- DNAString(as.character(aln['ASV_1']))
  ASV1_seq_no_gaps <- as.character(RemoveGaps(aln['ASV_1'],removeGaps = "all"))
  first10nt <- stringr::str_sub(ASV1_seq_no_gaps,1,10)
  last10nt <- stringr::str_sub(ASV1_seq_no_gaps,-10,-1)
  s <- start(matchPattern(first10nt, ASV1))
  e <-end(matchPattern(last10nt, ASV1))
  alnTrim <- DNAStringSet(aln,start=s,end=e)
  seqlengths = width(RemoveGaps(alnTrim,
                                removeGaps = "all",
                                processors = 1))
  alnTrimFilt <- alnTrim[seqlengths > 250*.95]
  alnTrimFilt <- RemoveGaps(alnTrimFilt,
                            removeGaps = "common")
  print(c(length(aln),'seqs in alignment',length(alnTrimFilt),'seqs in trimmed alignment'))
  Biostrings::writeXStringSet(alnTrimFilt,file=aln_out)
}

#ref taxa phylogeny
asv_ref <- c(bacteroidales_ref,DNAStringSet(asv_fasta_filt,start=2)) #read in nuc fasta with ref gyrb seqs from gtdbtk
Biostrings::writeXStringSet(asv_ref,file=file.path(processing_folder,'phylogeny',"ASVs_filtered_ref.fasta")) #write ASVs with ref fasta
asv_ref_aln <- AlignTranslation(asv_ref)   #align based on AA sequences
Biostrings::writeXStringSet(asv_ref_aln,file=file.path(processing_folder,'phylogeny',"ASVs_filtered_ref_aln.fasta"))  #write alignment
trimAlign(file.path(processing_folder,'phylogeny','ASVs_filtered_ref_aln.fasta'), #trim alignment to where ASV1 starts and ends
          file.path(processing_folder,'phylogeny',"ASVs_filtered_ref_aln_trim.fasta"))
system(paste0('./FastTree-2.1.9 -nt -gtr <  ',
              file.path(processing_folder,'phylogeny',"ASVs_filtered_ref_aln_trim.fasta"), #generate phylogeny using fasttree
              ' > ',
              file.path(processing_folder,'phylogeny',"ASVs_filtered_ref_aln_trim.tre")))
asv_ref_tree <- ape::read.tree(file.path(processing_folder,'phylogeny',"ASVs_filtered_ref_aln_trim.tre")) #read in phylogeny
asv_tips <- asv_ref_tree$tip.label[str_detect(asv_ref_tree$tip.label,'ASV')] #find ASVs seqs
asv_ref_tree <- keep.tip(asv_ref_tree,asv_tips) #remove ref seqs
outgroup_mrca <- findMRCA(asv_ref_tree,F082_ASVs) #find outgroup
asv_ref_tree_rooted <- reroot(asv_ref_tree,outgroup_mrca) #reroot tree based on outgroup
ape::write.tree(asv_ref_tree_rooted,"tmp.tre")
system(paste0("sed 's/Root;/;/g'  tmp.tre > ",
              file.path(processing_folder,'phylogeny',"ASVs_filtered_ref_aln_trim_rooted.tre")))
system("rm tmp.tre")

#any differences in tree tips?
setdiff(asv_ref_tree$tip.label,asv_tree$tip.label)
setdiff(asv_tree$tip.label,asv_ref_tree$tip.label)
setdiff(taxa_names(ps_Bacteroidales),asv_ref_tree$tip.label) #ASV removed by generating phylogeny
setdiff(taxa_names(ps_Bacteroidales),asv_tree$tip.label)

#prune physeq to taxa in phylogeny
ps_Bacteroidales_final <- prune_taxa(asv_tree$tip.label,ps_Bacteroidales)

#output physeq to inputs_folder
inputs_folder <- 'results/gyrb/inputs'
dir.create(inputs_folder, recursive=TRUE)
#write taxonomy table
tax_table_tab <- as.data.frame(tax_table(ps_Bacteroidales_final)) %>% 
  rownames_to_column(var='ASV')
write.table(tax_table_tab,file.path(inputs_folder,'physeq_Bacteroidales_taxonomy.txt'),quote=F,row.names=F,sep='\t')
#write metadata
metadata_tab <- as.data.frame(sample_data(ps_Bacteroidales_final)) 
metadata_tab$X.SampleID <- rownames(metadata_tab)
write.table(metadata_tab,file.path(inputs_folder,'physeq_metadata_passing_samples.txt'),quote=F,row.names=F,sep='\t')
#write otu/asv table
otu_table <- as.data.frame(otu_table(ps_Bacteroidales_final)) 
write.table(otu_table,file.path(inputs_folder,'physeq_Bacteroidales_asv_tab.txt'),quote=F,row.names=T,col.names=NA,sep='\t')
#write fasta
asv_fasta_final <- asv_fasta[names(asv_fasta) %in% taxa_names(ps_Bacteroidales_final)]
Biostrings::writeXStringSet(asv_fasta_final,file=file.path(inputs_folder,"physeq_Bacteroidales_asv.fasta")) 

#copy over phylogenies
system(paste('cp',file.path(processing_folder,'phylogeny',"ASVs_filtered_aln_trim_rooted.tre"),
              file.path(inputs_folder,'physeq_Bacteroidales_ASVs.tree'),sep=' '))
system(paste('cp',file.path(processing_folder,'phylogeny',"ASVs_filtered_ref_aln_trim_rooted.tre"),
             file.path(inputs_folder,'physeq_Bacteroidales_ASVs_ref.tree'),sep=' '))

#write moeller codiv seq files
system(paste0('cp results/gyrb/processing/moeller_sup/moeller_codiv* ',inputs_folder))
