#This script is a function called by DADA2_V4.R 
#performs DADA2 QC/trimming/error correction/chimera removal on single end reads
#inputs sample fastqs with primers trimmed for an individual V4 dataset
#outputs an ASVs tableand summary of reads removed at each filtering step.

DADA2_single <- function(forward_reads, sample.names, outdir, trunc_parameters){
#Defining output directory/fastqs
filtered_path <- paste0(outdir,"/filtered_reads" )
dir.create(filtered_path, recursive=TRUE)
filtered_forward <- file.path(filtered_path, paste0(sample.names, "_R1_filt.fastq"))
filtered_forward
#Quality trimming/filtering
filtered_out = filterAndTrim(fwd=forward_reads, filt=filtered_forward,
                             maxN=0, maxEE=2,
                             truncLen=trunc_parameters, 
                             compress=TRUE, verbose=TRUE,multithread =TRUE)

# selecting only samples that have reads after filtering
filtered_out_df = data.frame(filtered_out)
filtered_out_df$sample_names <- sample.names
sample.names <- filtered_out_df[,'sample_names'][filtered_out_df$reads.out>0]

filtered_forward <- file.path(filtered_path, paste0(sample.names, "_R1_filt.fastq"))

#Dereplication
derep_forward <- derepFastq(filtered_forward, verbose=TRUE)

#Learn error rates
tryCatch(
  expr = {
      errF <- readRDS(file.path(outdir,"error_rates_forward.rds"))
      message("Found error rate files")
    },
    error = function(e){
      message('Files not found, generating error rates')
      print(e)
      errF <- learnErrors(derep_forward,multithread=TRUE)
      saveRDS(errF,file.path(outdir,"error_rates_forward.rds"))
    })
errF <- readRDS(file.path(outdir,"error_rates_forward.rds"))

#Generate ASVs
dada_forward <- dada(derep_forward, err=errF, multithread=TRUE)
(names(dada_forward) <- sample.names)

#Construct sequence table
seqtab <- makeSequenceTable(dada_forward)
dim(seqtab)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print(sum(seqtab.nochim)/sum(seqtab))

#Save ASV table, change name to whats in metadata
saveRDS(seqtab.nochim, file.path(outdir,"seqtab.nochim.rds")) #Saves ASV table
seqtab.nochim <- readRDS(file.path(outdir,"seqtab.nochim.rds")) #Opens ASV table
table(nchar(colnames(seqtab.nochim)))

#Summary of reads removed at each filtering step
getN <- function(x) {sum(getUniques(x))}
denoisedF <- sapply(dada_forward, getN)
denoisedF <- data.frame(denoisedF) %>% rownames_to_column(var = 'sample_names')
nonchim <- rowSums(seqtab.nochim)
nonchim <- data.frame(nonchim) %>% rownames_to_column(var = 'sample_names')
track <- filtered_out_df %>% 
  full_join(denoisedF, on='sample_names') %>% 
  full_join(nonchim, on='sample_names') %>%
  select('sample_names','reads.in','reads.out','denoisedF','nonchim')
colnames(track) <- c('sample_names',"input", "filtered", "denoisedF", "nonchim")
write.table(track,file.path(outdir,"track.txt"),sep='\t',row.names = FALSE,quote=FALSE)
}