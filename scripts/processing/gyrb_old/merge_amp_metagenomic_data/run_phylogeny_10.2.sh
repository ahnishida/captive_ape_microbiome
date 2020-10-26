#!/usr/bin/env bash
for folder in GemRef_AlignMafft BtRef_AlignMafft noRef_AlignMafft;
do
cd /Volumes/AHN/captive_ape_microbiome/results/gyrb/processing/gyrb_amp_meta_datasets_10.2;
echo $folder;
transeq -sequence $folder/ref.fasta -outseq $folder/ref.faa;
cat $folder/ref.fasta $folder/ASVs_filtered.fasta > $folder/ASVs_filtered_ref.fasta;
cat $folder/ref.faa $folder/ASVs_filtered.faa> $folder/ASVs_filtered_ref.faa;
mafft $folder/ASVs_filtered_ref.faa > $folder/ASVs_filtered_ref_aln.faa;
tranalign -asequence  $folder/ASVs_filtered_ref.fasta -bsequence $folder/ASVs_filtered_ref_aln.faa -outseq $folder/ASVs_filtered_ref_aln.fasta;
done;