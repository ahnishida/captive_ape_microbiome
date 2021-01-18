# This script is a function called by DADA2_V4.R 
# performs DADA2 QC/trimming/error correction/chimera removal on paired end reads
# inputs sample fastqs with primers trimmed for an individual V4 dataset
# outputs an ASVs tableand summary of reads removed at each filtering step.

DADA2_paired <- function(forward_reads, reverse_reads, sample.names, outdir, trunc_parameters) {

#Defining output directory/fastqs
filtered_path <- paste0(outdir,"/filtered_reads" )
dir.create(filtered_path, recursive=TRUE)
filtered_forward <- file.path(filtered_path, paste0(sample.names, "_R1_filt.fastq"))
filtered_reverse <- file.path(filtered_path, paste0(sample.names, "_R2_filt.fastq"))

#Quality trimming/filtering
filtered_out = filterAndTrim(fwd=forward_reads, filt=filtered_forward, rev=reverse_reads, filt.rev=filtered_reverse,
                               truncLen=trunc_parameters, 
                               maxN=0, maxEE=2,
                               compress=TRUE, verbose=TRUE,multithread =TRUE)

#Dereplication
derep_forward <- derepFastq(filtered_forward, verbose=TRUE)
derep_reverse <- derepFastq(filtered_reverse, verbose=TRUE)

#Learn error rates
tryCatch(
  expr = {
      errF <- readRDS(file.path(outdir,"error_rates_forward.rds"))
      errR <- readRDS(file.path(outdir,"error_rates_reverse.rds"))
      message("Found error rate files")
    },
    error = function(e){
      message('Files not found, generating error rates')
      print(e)
      errF <- learnErrors(derep_forward,multithread=TRUE)
      saveRDS(errF,file.path(outdir,"error_rates_forward.rds"))
      errR <- learnErrors(derep_reverse,multithread=TRUE)
      saveRDS(errR,file.path(outdir,"error_rates_reverse.rds"))
    })

errF <- readRDS(file.path(outdir,"error_rates_forward.rds"))
errR <- readRDS(file.path(outdir,"error_rates_reverse.rds"))

#Generate ASVs
dada_forward <- dada(derep_forward, err=errF, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=errR, multithread=TRUE)
#Merge forward and reverse reads
merged <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, verbose=TRUE)
(names(merged) <- sample.names)
#Construct sequence table
seqtab <- makeSequenceTable(merged)
dim(seqtab)
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print(sum(seqtab.nochim)/sum(seqtab))

#Save ASV table, change name to whats in metadata
saveRDS(seqtab.nochim, file.path(outdir,"seqtab.nochim.rds")) #Saves ASV table
seqtab.nochim <- readRDS(file.path(outdir,"seqtab.nochim.rds")) #Opens ASV table
table(nchar(colnames(seqtab.nochim)))

#Summary of reads removed at each filtering step
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, sapply(dada_forward, getN),  rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track) %>% rownames_to_column(var='sample_name')
write.table(track,file.path(outdir,"track.txt"),sep='\t',row.names = FALSE,quote=FALSE) 
}


