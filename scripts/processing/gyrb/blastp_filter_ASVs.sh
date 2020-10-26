#!/usr/bin/env bash
input=$1
db=$2
output=$3
echo $output
blastp -query $input -db $db -outfmt "7 qseqid sacc pident qlen length evalue bitscore salltitles sseq qseq sstart send" -max_target_seqs 5 -out $output
	
