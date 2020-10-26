library(stringr)
library(phytools)

setwd('/Volumes/AHN/captive_ape_microbiome/results/gyrb/processing/gyrb_amp_meta_datasets_10.2')
TREE <- ape::read.tree("results_delta/ASVs_filtered_aln.tre")
outgroup_taxa = TREE$tip.label[str_detect(TREE$tip.label, 'p__Gemmatimonadota')]
mrca <- findMRCA(TREE,outgroup_taxa)
TREE_rooted <- reroot(TREE,mrca)
write.tree(TREE_rooted,"results_delta/ASVs_filtered_aln_rooted.tre")
