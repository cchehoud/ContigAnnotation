#!/bin/bash

#BLAST a protein query against all of the viral families, and append the name of each viral family to the last column of the output, creating only a single output file
quer=$1
out=$2 #Will append to provided output file
ref_db=$3
if [ "$quer" == "" ]; then echo "Must provide input file"; exit; fi
if [ "$out" == "" ]; then echo "Must provide output file"; exit; fi
if [ "$ref_db" == "" ]; then echo "Must provide reference viral families db files"; exit; fi

cp $quer temp.query.fasta
cat $ref_db/virus.families | while read fam
do
    blastp -query temp.query.fasta -db $ref_db/$fam -outfmt 7 -evalue 1e-10 -num_threads 20| grep -v "#" | awk -v fam="$fam" '{print $0 "\t" fam}' >> $out
done
rm temp.query.fasta
