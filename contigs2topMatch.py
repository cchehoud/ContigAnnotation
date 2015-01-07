#!/usr/bin/python

import argparse

def extract_match (filein):
    filein.seek(0)
    contigs2match = {}
    for line in filein:
        line = line.rstrip()
        values = line.split("\t")
        contig = values[0]
        match = values[1]
        if not (contig in contigs2match):
            contigs2match[contig] = match
    return(contigs2match)
            
if __name__ == '__main__':
    '''
    This program will take BLAST output file (format 6 with qseqid and stitle) and return a hash with the query and the top blast match.
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Top Match.')
    parser.add_argument('-i', '--input', dest='blastFile', required=True, type=file,
                         help='Required input: BLAST file in format 6 with qseqid and stitle.') 

    args = parser.parse_args()

    print "Extracting information from BLAST file:"
    contigs2count = extract_match(args.blastFile)

    print "ContigName\tTopMatch"
    for contig, match in contigs2count.items():
        print contig + "\t" + str(match)
