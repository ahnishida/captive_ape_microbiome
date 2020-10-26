cd /Volumes/AHN/captive_ape_microbiome

gyrb_ref_full_fasta=results/gyrb/processing/ref_gyrb_gtdbtk/gyrb_fastas/gtdbtk_gyrb_Bt.fasta
gyrb_ref_full_faa=results/gyrb/processing/ref_gyrb_gtdbtk/gyrb_fastas/gtdbtk_gyrb_Bt.faa

asv_fasta=results/gyrb/processing/gyrb_amp_meta_datasets_10.2/ASVs_filtered.fasta
asv_faa=results/gyrb/processing/gyrb_amp_meta_datasets_10.2/ASVs_filtered.faa

mkdir results/gyrb/processing/gyrb_amp_meta_datasets_10.2/phylogeny
asv_fasta_ref_full=results/gyrb/processing/gyrb_amp_meta_datasets_10.2/phylogeny/ASVs_filtered_alpha.fasta
asv_faa_ref_full=results/gyrb/processing/gyrb_amp_meta_datasets_10.2/phylogeny/ASVs_filtered_alpha.faa


transeq -frame 2 -sequence $asv_fasta -outseq $asv_faa
cat $gyrb_ref_full_fasta $asv_fasta > $asv_fasta_ref_full 
cat $gyrb_ref_full_faa $asv_faa > $asv_faa_ref_full 
mafft $asv_faa_ref_full > ${asv_faa_ref_full}.aln
tranalign -asequence  $asv_fasta_ref_full -bsequence ${asv_faa_ref_full}.aln -outseq ${asv_fasta_ref_full}.aln

ipython scripts/processing/gyrb_old/merge_amp_metagenomic_data/trim_alignment.ipynb ${asv_fasta_ref_full}.aln ${asv_fasta_ref_full}.aln.full 421 3117
ipython scripts/processing/gyrb_old/merge_amp_metagenomic_data/trim_alignment.ipynb ${asv_fasta_ref_full}.aln ${asv_fasta_ref_full}.aln.amp 727 1227

fasttree -nt -gtr <  ${asv_fasta_ref_full}.aln.amp >  ${asv_fasta_ref_full}.aln.amp.tree
fasttree -nt -gtr <  ${asv_fasta_ref_full}.aln.full >  ${asv_fasta_ref_full}.aln.full.tree

R -f scripts/processing/gyrb_old/merge_amp_metagenomic_data/root_alpha_tree.R
