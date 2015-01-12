#!/usr/bin/python

import argparse
import sys

def extract_name_length (filein):
    filein.seek(0)
    contigs2length = {}
    for line in filein:
        line = line.rstrip('\n')
        if line.startswith(">"):
            values = line.split(" ")
            contig_name = values[0].split(">")[1]
        else:
            length = len(line)
            contigs2length[contig_name] = length
    return(contigs2length)
    
if __name__ == '__main__':
    '''
    This program will take a FASTA file of contigs and extract the name and length. The format of the identifier is as follows ">ContigName ...". For example: ">contig-95_2 [followed by any or no text]". 
    '''

    parser = argparse.ArgumentParser(description='Contig Names to Length and Read Count.')
    parser.add_argument('-i', '--input', dest='contigFile', required=True, type=file,
                         help='Required input: contig file in fasta format.')

    args = parser.parse_args()

    contigs2length = extract_name_length(args.contigFile)
    print "ContigName\tLength"
    for contig, length in contigs2length.items():
        print contig + "\t" + str(length)
