#!/usr/bin/python

import argparse

def extract_name_length_readcount (filein):
    contigs2length = {}
    contigs2readcount = {}
    for line in filein:
        line = line.rstrip('\n')
        if line.startswith(">"):
            values = line.split(" ")
            contig_name = values[0].split(">")[1]
            length = values[1].split("length_")[1]
            read_count = values[2].split("read_count_")[1]
            contigs2length[contig_name] = length
            contigs2readcount[contig_name] = read_count
    return(contigs2length,contigs2readcount)
    
if __name__ == '__main__':
    '''
    This program will take a FASTA file of contigs (outputed by glimmer) and extract the name, length, and read counts. 
    '''

    parser = argparse.ArgumentParser(description='Contig Names to Length and Read Count.')
    parser.add_argument('-i', '--input', dest='contigFile', required=True, type=file,
                         help='Required input: contig file in fasta format.')

    args = parser.parse_args()

    print "Extracting information from contig file:"
    extract_name_length_readcount(args.contigFile)

