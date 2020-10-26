library(dada2); packageVersion("dada2");
library(phyloseq)
library(tidyverse)

setwd('/Volumes/AHN/captive_ape_microbiome')

metagenomic_path <- "Pasolli_2019/metagenomic_assemblies/DADA_input_files" 
filtered_path <- "Pasolli_2019/metagenomic_assemblies/DADA_output_files" 
dir.create(filtered_path, recursive=TRUE)

#Input files and sample names
(forward_reads <- sort(list.files(metagenomic_path, pattern=".fasta", full.names = TRUE)))
(forward_reads = forward_reads[1:10])

(sample.names <- sapply(strsplit(basename(forward_reads), ".fasta"), `[`, 1))

seqtab <- makeSequenceTable(forward_reads)


#Dereplication
derep_forward <- derepFastq(filtered_forward, verbose=TRUE)
#Learn error rates
errF <- learnErrors(derep_forward, multithread=FALSE)
saveRDS(errF,file.path(DADA2_path,"error_rates_forward.rds"))
#Generate ASVs
dada_forward <- dada(derep_forward, err=errF, multithread=FALSE)
(names(dada_forward) <- sample.names)
#Construct sequence table
seqtab <- makeSequenceTable(dada_forward)
dim(seqtab)
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)


#Save ASV table, change name to whats in metadata
metadata_file='data/gyrb.family_metagenomic/metadata/metadata.txt'
metadata <- read.csv(metadata_file,sep='\t')
as.character(metadata$Old_SampleID)

sample_names(otu_table(seqtab.nochim,taxa_are_rows=FALSE)) #verify that table imports in Phyloseq
getid <- function(old_id){
  new_id= metadata$X.SampleID[metadata$Old_SampleID==old_id]
  return(new_id)
}
(new_names<- sapply(strsplit(basename(row.names(seqtab.nochim)), "_"), `[`, 1))
(new_names <- as.character(unlist(lapply(new_names,getid))))
row.names(seqtab.nochim)<- new_names

#Save ASV table
saveRDS(seqtab.nochim, file.path(DADA2_path,"seqtab.nochim.rds")) #Saves ASV table
seqtab.nochim <- readRDS(file.path(DADA2_path,"seqtab.nochim.rds")) #Opens ASV table

#Summary of reads removed at each filtering step
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, sapply(dada_forward, getN),  rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track) %>% rownames_to_column(var='sample_name')
write.table(track,file.path(DADA2_path,"track.txt"),sep='\t',row.names = FALSE,quote=FALSE)

