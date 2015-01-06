#!/usr/bin/python

import sys
import argparse
import subprocess
import os
import circular
import contigs2length
import contigs2circular
import contigs2ORFs
import contig2ViralFamily
import contigs2count

def run_glimmer(contig_file):
    subprocess.check_call(["bash", "glimmer-wrapper.sh", contig_file.name])
 
def run_blast_viraldb(ORF_file, viral_Blast_file, ref_viralDB):
    subprocess.check_call(["bash", "blast.viral.families.sh", ORF_file, viral_Blast_file, ref_viralDB])

def run_blast_against_db(ORF_file, output_Blast_file, ref_blastDB):
    subprocess.check_call(["blastp", "-query", ORF_file, "-db", ref_blastDB, "-outfmt", "7", "-evalue", "1e-10", "-out", output_Blast_file])
     
def extract_annotations(contig_file_fh, circle_file_fh, orf_file_fh, blast_file_fh,integrase_Blast_fh,aclame_Blast_fh, vfdb_Blast_fh):
    (c2length, c2readcount) = contigs2length.extract_name_length_readcount(contig_file_fh)
    c_circular = contigs2circular.extract_circularity(circle_file_fh)
    c2ORFs = contigs2ORFs.extract_ORF_counts(orf_file_fh)
    c2viralORFs = contigs2count.extract_counts(blast_file_fh)
    c2familyName = contig2ViralFamily.extract_family_name(blast_file_fh)
    c2integrase = contigs2count.extract_counts(integrase_Blast_fh)
    c2aclame = contigs2count.extract_counts(aclame_Blast_fh)
    c2vfdb = contigs2count.extract_counts(vfdb_Blast_fh)
    return (annotation_table(c2length, c2readcount, c_circular,c2ORFs,c2viralORFs, c2familyName,c2integrase,c2aclame,c2vfdb))

def annotation_table(c2length, c2readcount, c_circular, c2ORFs, c2viralORFs, c2familyName,c2integrase,c2aclame,c2vfdb):
    table = []
    for contig, length in c2length.items():
        readcount = c2readcount[contig]
        if contig in c_circular:
            circular = "Yes"
        else:
            circular = "No"
        if contig in c2ORFs:
            nORFs = c2ORFs[contig]
        else:
            nORFs = 0
        if contig in c2viralORFs:
            nViralORFs = c2viralORFs[contig]
        else:
            nViralORFs = 0
        if contig in c2familyName:
            Family = c2familyName[contig]
        else:
            Family = "NA"
        if contig in c2integrase:    
            integrase = c2integrase[contig]
        else:
            integrase = 0
        if contig in c2aclame:
            aclame = c2aclame[contig]
        else:
            aclame = 0
        if contig in c2vfdb:
            vfdb = c2vfdb[contig]
        else:
            vfdb = 0

        table.append([contig, length, readcount, circular, nORFs, nViralORFs, Family, integrase,aclame, vfdb])
    return(table)

if __name__ == '__main__':
        '''
        This program will take a FASTA file of contigs in nucleotide form and output a table contianining circulairty and ORF information. The format of the table is as follows: contigName, length, readCount, circular, numOfORFS, numOfORFSMatchingViralFamily, bestViralFamilyClassification.
        '''

        parser = argparse.ArgumentParser(description='Contig Annotation.')
        parser.add_argument('-i', '--input', dest='contigFile', required=True, type=file,
                            help='Required input: contig file in fasta format.')
        parser.add_argument('-o', '--output_contig_annotation_table', dest='outputFile', required=False,
                            type=str, help='Optional output file name.')
#        parser.add_argument('-p' '--proteinDB', dest='ref_protein_db', required=True, type=file,
#                            help='Required input: comma delimited file of protein reference db name and path.')
#        parser.add_argument('-n' '--nucleotideDB', dest='ref_nucleotide_db', required=True, type=file,
#                            help='Required input: comma delimited file of nucleotide reference db name and path.')
        parser.add_argument('-r', '--ref_viral_protein_DB', dest='ref_viral', required=True,
                            type=str,
                            help='Required input: location of db of viral family proteins.')        
        parser.add_argument('-g', '--integrase', dest='ref_integraseDB', required=True, type=str,
                           help='Required input: location of db of integrase proteins.')
        parser.add_argument('-a', '--aclame', dest='ref_aclameDB', required=True, type=str,
                            help='Required input: location of db of ACLAME phage proteins.')
        parser.add_argument('-v', '--vfdb', dest='ref_virulencefactorDB', required=True, type=str,
                            help='Required input: location of db of Virulence Factors.')
        
        args = parser.parse_args()

        basename = os.path.splitext(args.contigFile.name)[0]
        
        print "Running Circularity check:"
        circle_file = basename + "_circularity.fasta"
        circle_fh = open(circle_file, 'w')
        kMerMin = 10
        kMerMax = 1000
        min_len = 3500
        circular.find_circular_by_kmer(args.contigFile, circle_fh, kMerMin, kMerMax, min_len)
        circle_fh.close()
        circle_fh = open(circle_file, 'r')
        
        print "Running Glimmer:" 
        run_glimmer(args.contigFile)

        ORF_file = basename + ".fastp"
        ORF_file_fh = open(ORF_file, 'r')
        viral_Blast_file = basename + "_viral_blastout.txt"
        print "Running BLAST on viral db:"
        run_blast_viraldb(ORF_file, viral_Blast_file, args.ref_viral)
        viral_Blast_fh = open(viral_Blast_file, 'r')

        ## set of filehandles
        #for line in args.ref_protein_db:
        #    values = line.split(",")
        #    refDB_name = values[0]
        #    refDB_path = values[1]
        #    output_file = basename + refDB_name
        #    print "Running BLAST on " + refDB_name + " db:"
        #    run_blast_against_db(ORF_file, output_file, refDB_path)
            
        integrase_Blast_file = basename + "_integrase_blastout.txt"
        print "Running BLAST on integrase db:"
        run_blast_against_db(ORF_file, integrase_Blast_file, args.ref_integraseDB)
        integrase_Blast_fh = open(integrase_Blast_file, 'r')
        
        aclame_Blast_file = basename + "_aclame_blastout.txt"
        print "Running BLAST on ACLAME db:"
        run_blast_against_db(ORF_file, aclame_Blast_file, args.ref_aclameDB)
        aclame_Blast_fh = open(aclame_Blast_file, 'r')

        vfdb_Blast_file = basename + "_vfdb_blastout.txt"
        print "Running BLAST on VFDB db:"
        run_blast_against_db(ORF_file, vfdb_Blast_file, args.ref_virulencefactorDB)
        vfdb_Blast_fh = open(vfdb_Blast_file, 'r')

        table = extract_annotations(args.contigFile, circle_fh, ORF_file_fh, viral_Blast_fh, integrase_Blast_fh, aclame_Blast_fh, vfdb_Blast_fh)

        if (args.outputFile): 
            output = open(args.outputFile, 'w')
        else:
            output = sys.stdout
            
        for line in table:
            output.write('\t'.join(map(str, line)))
            output.write('\n')
 
