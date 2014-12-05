#!/usr/bin/python

import argparse

def extract_circularity (filein):
    contigs2circularity = {}
    for line in filein:
        line = line.rstrip('\n')
        if line.startswith(">"):
            values = line.split(" ")
            contig_name = values[0].split(">")[1]
            contigs2circularity[contig_name] = 'circular'
    return(contigs2circularity)
    
if __name__ == '__main__':
    '''
    This program will the output of the circularity program (fasta file of circular contigs) and return a hash of the contigs and their circularity. 
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Circularity.')
    parser.add_argument('-i', '--input', dest='circularContigsFile', required=True, type=file,
                         help='Required input: file of circular contigs.')

    args = parser.parse_args()

    print "Extracting information from contig file:"
    extract_circularity(args.circularContigsFile)

