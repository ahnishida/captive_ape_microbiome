#!/usr/bin/env bash

# generate reference gyrb fasta from gtdbtk database
# #ipython scripts/processing/gyrb/create_Bacteroidetes_GTDBTK_ref.ipynb

# process moeller et al 2016 supplemental data
# #ipython scripts/processing/gyrb/moeller_sup_codiv_clades.ipynb

#process amplicon data
##./scripts/processing/gyrb/gyrb_cutadapt.sh
##Rscript scripts/processing/gyrb/DADA2_gyrb_Bt.R
  #calls DADA2_single.R

# process metagenomic data
# ipython scripts/processing/gyrb/filter_gyrb_seqs_from_metagenomic_samples.ipynb
# ipython scripts/processing/gyrb/generate_metadata_for_metagenomic_samples.ipynb

#merge amplicon and metagenomic
Rscript scripts/processing/gyrb/merge_amplicon_metagenomic_datasets.R
  #calls idTaxa.R to assign taxonomy
  #calls blastp_filter_ASVs.sh to execute blastp
