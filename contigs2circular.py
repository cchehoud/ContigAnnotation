#!/usr/bin/python

import argparse

def extract_circularity (file_fh):
    circularcontigs = []
    
    for line in file_fh:
        line = line.rstrip('\n')
        if line.startswith(">"):
            values = line.split(" ")
            contig_name = values[0].split(">")[1]
            circularcontigs.append(contig_name)
    return(set(circularcontigs))
    
if __name__ == '__main__':
    '''
    This program will the output of the circularity program (fasta file of circular contigs) and return a set of the circular contigs.
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Circularity.')
    parser.add_argument('-i', '--input', dest='circularContigsFile', required=True, type=file,
                         help='Required input: file of circular contigs.')

    args = parser.parse_args()

    circularcontigs = extract_circularity(args.circularContigsFile)
    
    print "circularContigNames"
    for contig in circularcontigs:
        print contig

