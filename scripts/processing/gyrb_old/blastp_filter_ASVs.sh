#!/usr/bin/env bash
input=$1
db=$2
output=$3
echo $output
blastp -query $input -db $db -outfmt '6 qseqid sseqid pident length qlen evalue'  -max_target_seqs 1 -out $output
	
