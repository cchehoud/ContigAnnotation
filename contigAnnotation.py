#!/usr/bin/python

import argparse
import subprocess
import os
import circular

def run_glimmer(contig_file):
    subprocess.check_call(["bash", "glimmer-wrapper.sh", contig_file.name])

 
def run_blast_viraldb(ORF_file, viral_Blast_file, ref_viralDB):
    subprocess.check_call(["bash", "blast.viral.families.sh", ORF_file, viral_Blast_file, ref_viralDB])
    
if __name__ == '__main__':
        '''
        This program will take a FASTA file of contigs in nucleotide form and output a table contianining circulairty and ORF information. The format of the table is as follows: contigName, length, readCount, circular, numOfORFS, numOfORFSMatchingViralFamily, bestViralFamilyClassification.
        '''

        parser = argparse.ArgumentParser(description='Contig Annotation.')
        parser.add_argument('-i', '--input', dest='contigFile', required=True, type=file,
                            help='Required input: contig file in fasta format.')
        parser.add_argument('-r', '--ref_viral_protein_DB', dest='ref_viral', required=True,
                            type=str,
                            help='Required input: location of db of viral family proteins.') 
        
        args = parser.parse_args()

        basename = os.path.splitext(args.contigFile.name)[0]
        
        print "Running Circularity check:"
        circle_file = basename + "_circularity.fasta"
        circle_fh = open(circle_file, 'w')
        kMerMin = 10
        kMerMax = 1000
        min_len = 3500
        circular.find_circular_by_kmer(args.contigFile, circle_fh, kMerMin, kMerMax, min_len)

        
        print "Running Glimmer:" 
        run_glimmer(args.contigFile)

        ORF_file = basename + ".fastp"
        viral_Blast_file = basename + "_viral_blastout.txt"
        
        print "Running BLAST on viral db:"
        run_blast_viraldb(ORF_file, viral_Blast_file, args.ref_viral)

        print "Done"
