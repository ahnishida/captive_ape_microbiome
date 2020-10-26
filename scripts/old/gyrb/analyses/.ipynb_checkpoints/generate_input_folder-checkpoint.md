---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```bash
cd /Volumes/AHN/captive_ape_microbiome/
cp metadata/metadata_gyrb_amp_meta_passing_samples.txt results/gyrb/inputs/metadata_gyrb_amp_meta_passing_samples.txt
cd results/gyrb/

cp processing/gyrb_amp_meta/ASVs_filtered.fasta inputs/ASVs_filtered.fasta
cp processing/gyrb_amp_meta/ASVs_filtered_counts.tsv inputs/ASVs_filtered_counts.tsv
cp processing/gyrb_amp_meta/assign_taxonomy/ASVs_taxonomy.txt inputs/ASVs_taxonomy.txt

cp processing/gyrb_amp_meta/phylogeny/ASVs_filtered_ref_full.fasta.aln.full.rooted.tree inputs/ASVs_filtered_ref_full.tree

cp processing/moeller_sup/moeller_codiv_Bacteroidaceae.fna inputs/moeller_codiv_Bacteroidaceae.fna
cp processing/moeller_sup/moeller_codiv_HRclades.txt inputs/moeller_codiv_HRclades.txt
cp processing/moeller_sup/moeller_codiv_lin_Bt1.tree inputs/moeller_codiv_lin_Bt1.tree
cp processing/moeller_sup/moeller_codiv_lin_Bt2.tree inputs/moeller_codiv_lin_Bt2.tree
cp processing/moeller_sup/moeller_codiv_lin_Bt3.tree inputs/moeller_codiv_lin_Bt3.tree

```

```python

```

```python

```
