#!/bin/bash

#use cutadapt to trim primers from datasets

#trial cutadapt for one sample of each dataset shows only raymann (FandR primers),projectchimps (FandR primers), vangay (FandR primer), and moeller chimp (R primer)
#raw data 
infile=data/16s/16s_clayton_captive/raw_fastq/CZ-MB26_S26_L001 #just showing that primer trimming was already done
cutadapt -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT  -o test1.fastq -p test2.fastq $infile"_R1_001.fastq.gz" $infile"_R2_001.fastq.gz" 

infile=data/16s/16s_raymann_captive/raw_fastq/Orang114 #F and R primers need to be trimmed
cutadapt -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT  -o test1.fastq -p test2.fastq $infile"_R1_001.fastq.gz" $infile"_R2_001.fastq.gz" 

infile=data/16s/16s_smits_2017/raw_fastq/SRR5760852  #just showing that primer trimming was already done
cutadapt -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT  -o test1.fastq -p test2.fastq $infile"_1.fastq.gz" $infile"_2.fastq.gz" 

infile=data/16s/16s_vangay_2018/raw_fastq/CS.001  #F and R primers need to be trimmed
cutadapt -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT  -o test1.fastq -p test2.fastq $infile".R1.fastq.gz" $infile".R2.fastq.gz" 

infile=data/16s/16s_nishida_projectchimps/raw_fastq/ProjectChimps-A1-V4_S18_L001  #F and R primers need to be trimmed
cutadapt -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT  -o test1.fastq -p test2.fastq $infile"_R1_001.fastq.gz" $infile"_R2_001.fastq.gz" 

infile=data/16s/16s_goodrich/raw_fastq/ERR1382855 #just showing that primer trimming was already done
cutadapt -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT  -o test1.fastq -p test2.fastq $infile"_1.fastq.gz" $infile"_2.fastq.gz" 

#demultiplexed output
infile=results/16s/processing/16s_nishida_captive/demultiplexed_fastq/cp.bon.COLZ.1.16s #just showing that primer trimming was already done
cutadapt -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT  -o test1.fastq -p test2.fastq $infile"-read-1.fastq.gz" $infile"-read-2.fastq.gz" 

cutadapt -g GTGCCAGCCGCCGCGGTAA -o test1.fastq results/16s/processing/16s_moeller_wild/gor_bon/demultiplexed_fastq/wd.bon.LA.a187.16s.fastq.gz #just showing that primer trimming was already done
cutadapt -a ATTAGAWACCCBNGTAGTCC -o test1.fastq results/16s/processing/16s_moeller_wild/gor_bon/demultiplexed_fastq/wd.bon.LA.a187.16s.fastq.gz

cutadapt -g GTGCCAGCCGCCGCGGTAA -o test1.fastq  results/16s/processing/16s_moeller_wild/chimp/demultiplexed_fastq/wd.chi.GM.160.16s-read-1.fastq.gz
cutadapt -a ATTAGAWACCCBNGTAGTCC -o test1.fastq results/16s/processing/16s_moeller_wild/chimp/demultiplexed_fastq/wd.chi.GM.160.16s-read-1.fastq.gz #R primers need to be trimmed

rm test1.fastq; rm test2.fastq

# cutadapt raymann, trim forward and reverse
mkdir -pv results/16s/processing/16s_raymann_captive/cutadapt_fastq/
for f in data/16s/16s_raymann_captive/raw_fastq/*_R1_001.fastq.gz
do
sample=${f#data/16s/16s_raymann_captive/raw_fastq/}
sample=${sample%_R1_001.fastq.gz}
cutadapt \
	-g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed \
	-o results/16s/processing/16s_raymann_captive/cutadapt_fastq/"$sample"_R1_001.fastq.gz \
	-p results/16s/processing/16s_raymann_captive/cutadapt_fastq/"$sample"_R2_001.fastq.gz \
	data/16s/16s_raymann_captive/raw_fastq/"$sample"_R1_001.fastq.gz  data/16s/16s_raymann_captive/raw_fastq/"$sample"_R2_001.fastq.gz > \
	results/16s/processing/16s_raymann_captive/cutadapt_fastq/"$sample"_report.txt
echo $sample
done

# cutadapt nishida project chimps, trim forward and reverse
mkdir -pv results/16s/processing/16s_nishida_projectchimps/cutadapt_fastq/
for f in data/16s/16s_nishida_projectchimps/raw_fastq/*_R1_001.fastq.gz
do
sample=${f#data/16s/16s_nishida_projectchimps/raw_fastq/}
sample=${sample%_R1_001.fastq.gz}
cutadapt \
	-g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed \
	-o results/16s/processing/16s_nishida_projectchimps/cutadapt_fastq/"$sample"_R1_001.fastq.gz \
	-p results/16s/processing/16s_nishida_projectchimps/cutadapt_fastq/"$sample"_R2_001.fastq.gz \
	data/16s/16s_nishida_projectchimps/raw_fastq/"$sample"_R1_001.fastq.gz  data/16s/16s_nishida_projectchimps/raw_fastq/"$sample"_R2_001.fastq.gz > \
	results/16s/processing/16s_nishida_projectchimps/cutadapt_fastq/"$sample"_report.txt
echo $sample
done

# cutadapt vangay, trim forward and reverse
mkdir -pv results/16s/processing/16s_vangay_2018/cutadapt_fastq/
for f in data/16s/16s_vangay_2018/raw_fastq/*.R1.fastq.gz
do
sample="${f#data/16s/16s_vangay_2018/raw_fastq/}"
sample=${sample%.R1.fastq.gz}
cutadapt \
    -g GTGCCAGCCGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed \
	-o results/16s/processing/16s_vangay_2018/cutadapt_fastq/"$sample".R1.fastq.gz \
	-p results/16s/processing/16s_vangay_2018/cutadapt_fastq/"$sample".R2.fastq.gz \
	data/16s/16s_vangay_2018/raw_fastq/"$sample".R1.fastq.gz  data/16s/16s_vangay_2018/raw_fastq/"$sample".R2.fastq.gz > \
	results/16s/processing/16s_vangay_2018/cutadapt_fastq/"$sample"_report.txt
echo $sample
done

#cutadapt moeller chimp, trim reverse only
mkdir -pv results/16s/processing/16s_moeller_wild/chimp/cutadapt_fastq/
for fastq in results/16s/processing/16s_moeller_wild/chimp/demultiplexed_fastq/*-read-1.fastq.gz
do
sample="${fastq#results/16s/processing/16s_moeller_wild/chimp/demultiplexed_fastq/}"
sample="${sample%-read-1.fastq.gz}"
cutadapt \
    -a ATTAGAWACCCBNGTAGTCC --overlap 20 --discard-untrimmed \
	-o results/16s/processing/16s_moeller_wild/chimp/cutadapt_fastq/"$sample"-read-1.fastq.gz \
	results/16s/processing/16s_moeller_wild/chimp/demultiplexed_fastq/"$sample"-read-1.fastq.gz > \
	results/16s/processing/16s_moeller_wild/chimp/cutadapt_fastq/"$sample"_report.txt
echo $sample
done

