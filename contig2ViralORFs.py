#!/usr/bin/python

import argparse

def extract_ORF_counts (filein):
    contigs2ORFcount = {}
    contig_orfnames = set()
    for line in filein:
        values = line.split("\t")
        contig_orf = values[0]
        if not (contig_orf in contig_orfnames):
            contig_orfnames.add(contig_orf)
            contig_name = contig_orf.split(".")[0]
            orf_name = contig_orf.split(".")[1]
            if contig_name in contigs2ORFcount:
                contigs2ORFcount[contig_name] += 1
            else:
                contigs2ORFcount[contig_name] = 1
    return(contigs2ORFcount)
            
if __name__ == '__main__':
    '''
    This program will take BLAST output file (format 7 without comments), as outputed by blast.viral.families.sh and extract the number of viral ORFs.
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Number of Viral ORFS.')
    parser.add_argument('-i', '--input', dest='blastFile', required=True, type=file,
                         help='Required input: BLAST file in format 7 with no # lines.')

    args = parser.parse_args()

    print "Extracting information from BLAST file:"
    contigs2ORFcount = extract_ORF_counts(args.blastFile)

    print "ContigName\tNumberofViralORFS"
    for contig, nORFs in contigs2ORFcount.items():
        print contig + "\t" + str(nORFs)
