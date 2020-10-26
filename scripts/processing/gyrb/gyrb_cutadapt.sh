#uses cutadapt to trim forward primer from each of the samples, from each of the three datasets
#input: raw fastqs in data folder
#output: trimmed fastqs in results folder

folder=gyrb_nishida_projectchimps
mkdir -pv results/gyrb/processing/$folder/cutadapt_fastq/
for f in data/gyrb/$folder/raw_fastq/*_R1_001.fastq.gz
do
sample=${f#data/gyrb/$folder/raw_fastq/}
sample=${sample%_R1_001.fastq.gz}
cutadapt \
	-g CGGAGGTAARTTCGAYAAAGG  --overlap 21 -e .15 --discard-untrimmed \
	-o results/gyrb/processing/$folder/cutadapt_fastq/"$sample"_R1_001.fastq.gz \
	data/gyrb/$folder/raw_fastq/"$sample"_R1_001.fastq.gz >  \
	results/gyrb/processing/$folder/cutadapt_fastq/"$sample"_R1_001_report.txt
echo $sample
done

folder=gyrb_nishida_captive_wild
mkdir -pv results/gyrb/processing/$folder/cutadapt_fastq/
for f in data/gyrb/$folder/raw_fastq/*_R1_001.fastq.gz
do
sample=${f#data/gyrb/$folder/raw_fastq/}
sample=${sample%_R1_001.fastq.gz}
cutadapt \
	-g CGGAGGTAARTTCGAYAAAGG  --overlap 21 -e .15 --discard-untrimmed \
	-o results/gyrb/processing/$folder/cutadapt_fastq/"$sample"_R1_001.fastq.gz \
	data/gyrb/$folder/raw_fastq/"$sample"_R1_001.fastq.gz >  \
	results/gyrb/processing/$folder/cutadapt_fastq/"$sample"_R1_001_report.txt
echo $sample
done

folder=gyrb_moeller_wild
mkdir -pv results/gyrb/processing/$folder/cutadapt_fastq/
for f in data/gyrb/$folder/raw_fastq/*_R1_001.fastq.gz
do
sample=${f#data/gyrb/$folder/raw_fastq/}
sample=${sample%_R1_001.fastq.gz}
cutadapt \
	-g CGGAGGTAARTTCGAYAAAGG  --overlap 21 -e .15 --discard-untrimmed \
	-o results/gyrb/processing/$folder/cutadapt_fastq/"$sample"_R1_001.fastq.gz \
	data/gyrb/$folder/raw_fastq/"$sample"_R1_001.fastq.gz >  \
	results/gyrb/processing/$folder/cutadapt_fastq/"$sample"_R1_001_report.txt
echo $sample
done