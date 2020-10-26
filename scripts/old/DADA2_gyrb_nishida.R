library(dada2); packageVersion("dada2");
library(phyloseq)
library(tidyverse)

setwd('/Volumes/AHN/captive_ape_microbiome')

rawdata_path <- "results/gyrb.family_nishida_captive_wild_human/cutadapt_fastq/Bif" 
DADA2_path <- "results/gyrb.family_nishida_captive_wild_human/DADA2/Bif" 
filtered_path <- "results/gyrb.family_nishida_captive_wild_human/DADA2/Bif/filtered_reads" 
dir.create(filtered_path, recursive=TRUE)

#Input files and sample names
(forward_reads <- sort(list.files(rawdata_path, pattern=".fastq", full.names = TRUE)))
#(forward_reads = forward_reads[1:10])
(sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1))

#View read quality
plotQualityProfile(forward_reads[1:2])

#Quality trimming/filtering
filtered_forward <- file.path(filtered_path, paste0(sample.names, "_filtered.fastq"))
filtered_out = filterAndTrim(fwd=forward_reads, filt=filtered_forward, 
                             truncLen=c(250), 
                             maxN=0, maxEE=2,
                             compress=TRUE, verbose=TRUE)

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

#Save ASV table
sample_names(otu_table(seqtab.nochim,taxa_are_rows=FALSE)) #verify to see table imports in phyloseq
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


rawdata_path <- "results/gyrb.family_nishida_captive_wild_human/cutadapt_fastq/Bt" 
DADA2_path <- "results/gyrb.family_nishida_captive_wild_human/DADA2/Bt" 
filtered_path <- "results/gyrb.family_nishida_captive_wild_human/DADA2/Bt/filtered_reads" 
dir.create(filtered_path, recursive=TRUE)

#Input files and sample names
(forward_reads <- sort(list.files(rawdata_path, pattern=".fastq", full.names = TRUE)))
#(forward_reads = forward_reads[1:10])
(sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1))

#View read quality
plotQualityProfile(forward_reads[1:2])

#Quality trimming/filtering
filtered_forward <- file.path(filtered_path, paste0(sample.names, "_filtered.fastq"))
filtered_out = filterAndTrim(fwd=forward_reads, filt=filtered_forward, 
                             truncLen=c(250), 
                             maxN=0, maxEE=2,
                             compress=TRUE, verbose=TRUE)

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

#Save ASV table
sample_names(otu_table(seqtab.nochim,taxa_are_rows=FALSE)) #verify to see table imports in phyloseq
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


