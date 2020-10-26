library(dada2); packageVersion("dada2")
library(phyloseq)
library(DECIPHER)
library(tidyverse)
library(phytools)
library(seqinr)
setwd('/Volumes/AHN/captive_ape_microbiome')
gyrb_ref_faa = 'ref_seqs/gtdbtk_gyrb.faa'
gyrb_ref_fasta = 'ref_seqs/gtdbtk_gyrb.fasta'
gyrb_ref <- readDNAStringSet(gyrb_ref_fasta,format = 'fasta')
s <- strsplit(names(gyrb_ref), ";")
phylum <- sapply(s, `[`, 2)
class <- sapply(s, `[`, 3)
order <- sapply(s, `[`, 4)
genus <- sapply(s, `[`, 6)
table(class[phylum == "p__Bacteroidota"]) 
table(order[class == "c__Bacteroidia"]) 

names(gyrb_ref) <- gsub(';','_',names(gyrb_ref))
names(gyrb_ref) <- gsub(' ','__',names(gyrb_ref))

Bt <- gyrb_ref[order %in% c('o__Bacteroidales')]
outgroup <- gyrb_ref[order %in% c('o__AKYH767','o__NS11-12g')]
Biostrings::writeXStringSet(outgroup,'ref_seqs/Bt_outgroup_only.fasta')

Bt_out <- c(Bt,outgroup)

Biostrings::writeXStringSet(Bt_out,'ref_seqs/Bt_out.fasta')
Bt_out_aln <- AlignTranslation(Bt_out)

Biostrings::writeXStringSet(Bt_out_aln,'ref_seqs/Bt_out_aln.fasta')
system(paste0('./FastTree-2.1.9 -nt -gtr <  ref_seqs/Bt_out_aln.fasta > ref_seqs/Bt_out_aln.tree'))
TREE <- ape::read.tree("ref_seqs/Bt_out_aln.tree")
outgroup_taxa = TREE$tip.label[str_detect(TREE$tip.label, c('o__AKYH767','o__NS11-12g'))]
mrca <- findMRCA(TREE,outgroup_taxa)
TREE_rooted <- reroot(TREE,mrca)
write.tree(TREE_rooted,"ref_seqs/Bt_out_aln_rooted.tree")
system("sed 's/Root;/;/g' ref_seqs/Bt_out_aln_rooted.tree > ref_seqs/Bt_out_aln_rooted2.tree")
system("mv ref_seqs/Bt_out_aln_rooted2.tree ref_seqs/Bt_out_aln_rooted.tree")

#compare to other ref
