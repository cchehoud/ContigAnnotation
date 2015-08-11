#!/usr/bin/python

import argparse

def extract_counts (filein):
    if filein is None:
        return {}
    filein.seek(0)
    contigs2count = {}
    contig_orfnames = set()
    for line in filein:
        line = line.rstrip()
        contig_orf = line.split("\t")[0]
        if not (contig_orf in contig_orfnames):
            contig_orfnames.add(contig_orf)
            contig_name = contig_orf.split(".")[0]
            if contig_name in contigs2count:
                contigs2count[contig_name] += 1
            else:
                contigs2count[contig_name] = 1
    return(contigs2count)
            
if __name__ == '__main__':
    '''
    This program will take BLAST output file (format 7 with or without comments) and return a hash with the counts and the number of matches. 
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Number of Matches.')
    parser.add_argument('-i', '--input', dest='blastFile', required=True, type=file,
                         help='Required input: BLAST file in format 7 with or without # lines.')

    args = parser.parse_args()

    print "Extracting information from BLAST file:"
    contigs2count = extract_counts(args.blastFile)

    print "ContigName\tCounts"
    for contig, count in contigs2count.items():
        print contig + "\t" + str(count)
