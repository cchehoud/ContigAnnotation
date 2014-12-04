#!/bin/bash

seqs=$1

basename=${1%.fasta}
basename=${basename%.fa}
basename=${basename%.scafSeq}
echo $basename

cat $seqs | while read line; do if [ "${line:0:1}" == ">" ]; then echo ${line%%' '*}; else echo $line; fi; done > $basename.temp

echo "Finding orfs from $seqs and writing to $basename.orfs"
cat $basename.temp | tigr-glimmer build-icm $basename.icm
tigr-glimmer glimmer3 -g 100 $basename.temp $basename.icm $basename
cat $basename.predict | while read line
	do if [ "${line:0:1}" == ">" ]
		then seqname=${line#'>'}
		else
		orf="$seqname.${line%%' '*}"
		coords="${line#*' '}"
		echo -e "$orf\t$seqname\t$coords"
		fi
	done > $basename.predict.formatted
#cat $basename.predict.formatted
tigr-glimmer multi-extract -l 100 --nostop $basename.temp $basename.predict.formatted > $basename.genes
echo "source('translateFasta.R'); translateFasta('$basename.genes','$basename.fastp')"|R --no-save --no-restore --silent
