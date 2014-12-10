#!/usr/bin/python

import argparse
import sys

def extract_name_length_readcount (filein):
    filein.seek(0)
    contigs2length = {}
    contigs2readcount = {}
    for line in filein:
        line = line.rstrip('\n')
        if line.startswith(">"):
            values = line.split(" ")
            if not len(values) == 3:
                    print("ERROR: Contig file is not in the correct format.")
                    sys.exit(1)
            contig_name = values[0].split(">")[1]
            length = values[1].split("length_")[1]
            read_count = values[2].split("read_count_")[1]
            contigs2length[contig_name] = length
            contigs2readcount[contig_name] = read_count
    return(contigs2length,contigs2readcount)
    
if __name__ == '__main__':
    '''
    This program will take a FASTA file of contigs (outputed by IDBA-UD) and extract the name, length, and read counts. The format of the identifier is as follows ">ContigName length read_count". For example: ">contig-95_2 length_33501 read_count_391116".
    '''

    parser = argparse.ArgumentParser(description='Contig Names to Length and Read Count.')
    parser.add_argument('-i', '--input', dest='contigFile', required=True, type=file,
                         help='Required input: contig file in fasta format.')

    args = parser.parse_args()

    (contigs2length,contigs2readcount) = extract_name_length_readcount(args.contigFile)
    print "ContigName\tLength\tReadCount"
    for contig, length in contigs2length.items():
        readcount = contigs2readcount[contig]
        print contig + "\t" + length + "\t" + readcount
        
