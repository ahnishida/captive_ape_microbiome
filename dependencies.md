# Dependencies 
##### Re-pair reads in disordered fastq (repair.sh from bbmap v38.70)
##### Demultiplex fastqs with barcode-splitter(v0.18.6) or demultiplex(v1.0.1)
##### Remove primers with cutadapt (v2.5)
##### Read trimming, error correction, chimera removal and generate ASV table with DADA2(v1.16.0)
##### Align sequences with mafft(v7.309)
##### Generate ASV phylogenies with FastTree(v2.1.9)
##### General data processing in R with tidyverse(v1.3.0), stringr(v1.4.0), reshape2(v1.4.4), rstatix(v0.6.0)
##### Processing microbiome datasets in R with phyloseq(v1.32.0), DECIPHER(v2.16.1), genefilter(v1.70.0), seqinr(v3.6-1), ape(v5.4-1), phytools(v0.7-70), zoo(v1.8-8), picante(v1.8.2)
##### Sequence analysis with prodigal(v2.6.3), hmmer(v3.3), transeq and transalign from EMBOSS(v6.6.0.0), BLAST(v2.9.0)
##### General data processing in python with pandas, numpy
##### Sequence analysis and defining HR clades in python with Biopython(v1.77), ete3(v3.1.1)
##### Statistical analysis in R with PMCMRplus(1.5.1), RVAideMemoire(v0.9-78), broom(0.7.1)
##### Data visualization in R with ggplot2(v3.3.3), ggtree(v2.2.4), cowplot(v1.1.0), RColorBrewer(v1.1-2)

### Dependencies for analyses-only
#### functions.ipynb, 16S_ASV_sharing_HRtype.ipynb, gyrb_hr_clades.ipynb
Python modules: ete3(v3.1.1),pandas,numpy

#### gyrb_visualize_all.R, 16s_analyses_all.R
Rlibraries: ape(v5.4-1), ggplot2(v3.3.3), ggtree(v2.2.4) ,tidyverse(v1.3.0), phytools(v0.7-70), cowplot(v1.1.0)
phyloseq(v1.32.0), tidyverse(v1.3.0),picante(v1.8.2), PMCMRplus(1.5.1), cowplot(v1.1.0),
ggplot2(v3.3.3), reshape2(v1.4.4), RVAideMemoire(v0.9-78), RColorBrewer(v1.1-2), rstatix(v0.6.0), broom(0.7.1)

### Dependencies for processing
#### V4_processing
##### V4_demultiplex.sh
demultiplex(v1.0.1), repair.sh from bbmap (v38.70), barcode-splitter(v0.18.6)

##### V4_cutadapt.sh 
cutadapt(v2.5)

##### DADA2_V4.R, DADA2_single.R, DADA2_paired.R 
DADA2(v1.16.0), phyloseq(v1.32.0), tidyverse(v1.3.0), 

##### merge_DADA_to_Phyloseq_V4.R
Rpackages: DADA2(v1.16.0), phyloseq(v1.32.0), genefilter(v1.70.0), tidyverse(v1.3.0)
seqinr(v3.6-1), ape(v5.4-1), phytools(v0.7-70)
programs: FastTree(v2.1.9), mafft(v7.309)

#### gyrb_processing
##### create_Bacteroidetes_GTDBTK_ref.ipynb: 
programs=prodigal(v2.6.3),hmmer(v3.3),transeq from EMBOSS(v6.6.0.0),BLAST(v2.9.0)
python packages: Biopython(v1.77)

##### moeller_sup_codiv_clades.ipynb
python packages: Biopython(v1.77), ete3(v3.1.1),
programs: transeq and transalign from EMBOSS(v6.6.0.0),cutadapt(v2.5), FastTree(v2.1.9), mafft(v7.309)

##### gyrb_cutadapt.sh
cutadapt(v2.5)

##### DADA2_gyrb_Bt.R,DADA2_single.R
DADA2(v1.16.0), phyloseq(v1.32.0), tidyverse(v1.3.0), 

##### filter_gyrb_seqs_from_metagenomic_samples.ipynb
Biopython(v1.77), BLAST(v2.9.0)

##### merge_amplicon_metagenomic_datasets.R
Rpackages: DADA2(v1.16.0), phyloseq(v1.32.0), tidyverse(v1.3.0), DECIPHER(v2.16.1),genefilter(v1.70.0) ,seqinr(v3.6-1), phytools(v0.7-70), stringr(v1.4.0), zoo(v1.8-8) 
programs: BLAST(v2.9.0), FastTree(v2.1.9)
scripts: blastp_filter_ASVs.sh, idTaxa.R


