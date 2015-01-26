#!/usr/bin/python

import argparse
import sys

from Bio import SeqIO

def extract_name_length(filein):
    """ length on seq.

        creates a map: key is read/seq id, value is number of nucleotides(length)

        Raises:
            ValueError if ID has been already seen
    """
    contigs2length = {}
    filein.seek(0)
    for record in SeqIO.parse(filein, "fasta"):
        contig_name = record.id
        if contig_name in contigs2length:
            raise ValueError("FASTA file has 2 identical ids: " + contig_name)
        length = len(record)
        contigs2length[contig_name] = length
    return contigs2length
    
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
