#!/usr/bin/env bash
#./scripts/processing/16s/V4_demultiplex.sh
#./scripts/processing/16s/V4_cutadapt.sh
#Rscript scripts/processing/16s/DADA2_V4.R
#calls DADA2_single or paired depending on dataset
Rscript scripts/processing/16s/merge_DADA_to_Phyloseq_V4.R
