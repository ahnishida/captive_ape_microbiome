library(seqinr)
library(tidyverse)
setwd('/Volumes/AHN/captive_ape_microbiome')

#assign taxonomy using GTDBTK database
#fileppaths to gyrb ref db
prot_db='results/gyrb/processing/ref_gyrb_gtdbtk/blastdb/gtdbtk_gyrb.faa'
nucl_db='results/gyrb/processing/ref_gyrb_gtdbtk/blastdb/gtdbtk_gyrb.fasta'


assign_taxonomy <- function(asv_fasta_file,asv_faa_file,outfolder,pident_cutoff,length_cutoff) {
  
#create outfolder
dir.create(outfolder, recursive=TRUE)

#run blastn and blastp
system(paste0('blastn -query ',asv_fasta_file,
              ' -db ',nucl_db,
              ' -outfmt "7 qseqid sacc pident qlen length evalue bitscore salltitles sseq" ',
              ' -max_target_seqs 1 ',
              ' -out ',outfolder,"/ASVs_blastn.txt"))
system(paste0('blastp -query ',asv_faa_file,
              ' -db ',prot_db,
              ' -outfmt "7 qseqid sacc pident qlen length evalue bitscore salltitles sseq qseq sstart send" ',
              ' -max_target_seqs 1 ',
              ' -out ',outfolder,"/ASVs_blastp.txt"))

#start tax table with column of all ASVs, 
#(some will return any hits throught blast search)
asv_fasta <- read.fasta(asv_fasta_file,as.string = TRUE)
ASV <- names(asv_fasta)
ASV_tax_table <- data.frame(ASV)

#identity asvs that when translated have stop codons
asv_faa <- read.fasta(asv_faa_file,as.string = TRUE,seqonly=TRUE)
(ASVs_w_stop_codons <- ASV[str_detect(asv_faa,'\\*')])

#read in blastp results and filter
blastp <- read.table(file.path(outfolder,"ASVs_blastp.txt"))
colnames(blastp) <- c('ASV','sacc','pident','qlen','length','evalue','bitscore','genome','salltitles','species_id','sseq','qseq','sstart','send')
blastp$ASV <- str_sub(blastp$ASV, 1, str_length(blastp$ASV)-2)
blastp <- blastp %>% separate(salltitles,sep=';',into=c('Domain','Phylum','Class','Order',
                                                        'Family','Genus','Species')) %>%
  filter(pident>pident_cutoff & length>(83*length_cutoff))  %>%
  filter(!ASV %in% ASVs_w_stop_codons) %>%
  filter(Order=='o__Bacteroidales') %>%
  filter(sstart > 110 & sstart < 120) %>%
  filter(send > 185 & send < 202) %>%
  select('ASV','Phylum','Class','Order','Family')

#read in blastp results and filter
blastn <- read.table(file.path(outfolder,"ASVs_blastn.txt"))
colnames(blastn) <- c('ASV','sacc','pident','qlen','length','evalue','bitscore','genome','salltitles','species_id','sseq')
blastn <- blastn %>% separate(salltitles,sep=';',into=c('Domain','Phylum','Class','Order',
                                                        'Family','Genus','Species')) %>% 
  filter(pident>pident_cutoff & length>(250*length_cutoff)) %>% 
  mutate(Genus_pident_cutoff = ifelse(pident>80,Genus,'Unassigned')) %>% 
  mutate(Species_pident_cutoff = ifelse(pident>95,Species,'Unassigned')) %>%
  select('ASV','Genus','Species') 

tax_table <- ASV_tax_table %>% 
  left_join(blastp, by='ASV') %>% 
  left_join(blastn, by='ASV') %>%
  mutate_each(funs(replace(., which(is.na(.)), 'Unassigned')))

#removes ASVs that didn't return a blastp hit passing filters
tax_table_filt <- tax_table %>% filter(Family != 'Unassigned') 
write.table(tax_table_filt, file=file.path(outfolder,"ASVs_taxonomy.txt"),sep='\t',quote=F,row.names=F)

}

#inputs test
asv_fasta_file = 'results/gyrb/processing/gyrb_amp_datasets/ASVs_all.fasta'
asv_faa_file = 'results/gyrb/processing/gyrb_amp_datasets/ASVs_all.faa'
outfolder = 'results/gyrb/processing/gyrb_amp_datasets/assign_taxonomy'
pident_cutoff <- 50
length_cutoff <- .8
#assign_taxonomy(asv_fasta_file,asv_faa_file,outfolder,pident_cutoff,length_cutoff) 
  