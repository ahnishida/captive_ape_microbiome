#!/usr/bin/env bash
#gunzip results/gyrb/inputs/physeq_Bacteroidales_asv_tab.txt.gz
#set current working dir in scripts 
ipython scripts/analyses/16S_ASV_sharing_HRtype.ipynb
  #loads functions.ipynb to assign HR ASVs
Rscript scripts/analyses/16S_analyses_all.R
ipython scripts/analyses/gyrb_hr_clades.ipynb
Rscript scripts/analyses/gyrb_visualize_all.R
