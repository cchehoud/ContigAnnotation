#!/usr/bin/python

import argparse

def extract_ORF_counts (filein):
    if filein is None:
        return {}
    filein.seek(0)
    contigs2ORFcount = {}
    for line in filein:
        line = line.rstrip('\n')
        if line.startswith(">"):
            values = line.split(".")
            contig_name = values[0].split(">")[1]
            if contig_name in contigs2ORFcount:
                contigs2ORFcount[contig_name] += 1
            else:
                contigs2ORFcount[contig_name] = 1
    return(contigs2ORFcount)
    
if __name__ == '__main__':
    '''
    This program will take a protein FASTA file of contigs (outputed by glimmer) and extract the number of ORFs per contig.
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Number of ORFS.')
    parser.add_argument('-i', '--input', dest='contigFile', required=True, type=file,
                         help='Required input: contig file in protein fasta format.')

    args = parser.parse_args()

    print "Extracting information from contig file:"
    contigs2ORFcount = extract_ORF_counts(args.contigFile)

    print "ContigName\tNumberofORFS"
    for contig, nORFs in contigs2ORFcount.items():
        print contig + "\t" + str(nORFs)
