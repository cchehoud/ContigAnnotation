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
import contigs2topMatch
import contigs2CDD

def run_CDD(basename, ORF_file, ref_CDD_DB):
    CDD_xml_file = basename + "_cdd.xml"
    CDD_out_file = basename + "_cdd.txt"
    print "Running CDD search:"
    command = "rpsblast -query " + ORF_file + " -evalue 0.01 -seg no -outfmt 5 -db " + ref_CDD_DB + " -num_threads 4 -out " + CDD_xml_file
    command2 = "rpsbproc -i " + CDD_xml_file + " -o " + CDD_out_file
    subprocess.check_call(command, shell=True)
    subprocess.check_call(command2, shell=True)
    return(open(CDD_out_file, 'r'))
        
def run_glimmer(contig_file):
    print "Running Glimmer:"
    subprocess.check_call(["bash", "glimmer-wrapper.sh", contig_file.name])
 
def run_blast_viraldb(basename, ORF_file, ref_viralDB):
    viralp_Blast_file = basename + "_viralp_blastout.txt"
    print "Running BLAST on viral proteins db:"
    subprocess.check_call(["bash", "blast.viral.families.sh", ORF_file, viralp_Blast_file, ref_viralDB])
    return(open(viralp_Blast_file, 'r'))
    
def run_blastp_against_db(basename, name, ORF_file, path):
    output_Blast_file = basename + "_" + name + "_blastout.txt"
    print "Running BLAST on " + name + " db:"
    command = "blastp -query " + ORF_file + " -db " + path + " -outfmt '6 qseqid' -evalue 1e-10 -num_threads 4 -max_target_seqs 1 -out " + output_Blast_file
    subprocess.check_call(command, shell=True)
    return(open(output_Blast_file, 'r'))
    
def run_blastn_against_db(basename, name, contig_file, path):
    output_Blast_file = basename + "_" + name + "_blastout.txt"
    print "Running BLAST on " + name + " db:"
    command = "blastn -query " + contig_file + " -db " + path + " -outfmt '6 qseqid stitle' -num_threads 4 -evalue 1e-10 -max_target_seqs 1 -out " + output_Blast_file
    subprocess.check_call(command, shell=True)
    return(open(output_Blast_file, 'r'))
         
def extract_annotations(contig_file_fh, circle_file_fh, orf_file_fh, viralp_Blast_fh, protein2fh, nucleotide2fh, cdd_fh): 
    c2length=contigs2length.extract_name_length(contig_file_fh)
    c_circular = contigs2circular.extract_circularity(circle_file_fh)
    c2ORFs = contigs2ORFs.extract_ORF_counts(orf_file_fh)
    c2viralORFs = contigs2count.extract_counts(viralp_Blast_fh)
    c2familyName = contig2ViralFamily.extract_family_name(viralp_Blast_fh)
    c2domains = contigs2CDD.extract_domain_counts(cdd_fh)
    
    proteinDBmatches = {}
    for proteinDB_name, fh in  protein2fh.items():
        proteinDBmatches[proteinDB_name] = contigs2count.extract_counts(fh) 

    nucleotideDBmatches = {}
    for nucleotideDB_name, fh in nucleotide2fh.items():
        nucleotideDBmatches[nucleotideDB_name] = contigs2topMatch.extract_match(fh)
    return (annotation_table(c2length, c_circular, c2ORFs, c2viralORFs, c2familyName, c2domains, proteinDBmatches, nucleotideDBmatches)) 

def annotation_table(c2length, c_circular, c2ORFs, c2viralORFs, c2familyName, c2domains, proteinDBmatches, nucleotideDBmatches): 
    table = []
    header = ["contig_name", "length", "circular", "nORFs", "nViralORFs", "family", "nDomains"]
    header += proteinDBmatches.keys() + nucleotideDBmatches.keys()
    table.append(header)
    for contig, length in c2length.items():
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
        if contig in c2domains:
            num_domains = c2domains[contig]
        else:
            num_domains = 0

        protien_list = extract_blast_hits(contig, proteinDBmatches, 0)
        nuc_list = extract_blast_hits(contig, nucleotideDBmatches, "NA") 
        row = [contig, length, circular, nORFs, nViralORFs, Family, num_domains]
        row += protien_list + nuc_list
        table.append(row)
    return(table)

def extract_blast_hits(contig_name, db_map, default_value):
    contig_value =[]
    # value here can be string (eg top match for contig) or integer (eg number of matches)
    for name_DB, contig_name2value in db_map.items(): 
        if contig_name in contig_name2value:
            match = contig_name2value[contig_name]
        else:
            match = default_value
        contig_value.append(match)
    return(contig_value)
        
def parse_arguments():
    parser = argparse.ArgumentParser(description='Contig Annotation.')
    parser.add_argument('-i', '--input', dest='contigFile', required=True, type=file,
                        help='Required input: contig file in fasta format.')
    parser.add_argument('-o', '--output_contig_annotation_table', dest='outputFile', required=False,
                        type=str, help='Optional output file name.')
    parser.add_argument('-d' '--databases', dest='databases', required=True, type=file,
                        help='Required input: configuration file(INI) with reference names and paths.')
    return(parser.parse_args())

def run_circular(basename, contig_file):
    print "Running Circularity check:"
    circle_file = basename + "_circularity.fasta"
    circle_fh = open(circle_file, 'w')
    kMerMin = 10
    kMerMax = 1000
    min_len = 3500
    circular.find_circular_by_kmer(contig_file, circle_fh, kMerMin, kMerMax, min_len)
    circle_fh.close()
    circle_fh = open(circle_file, 'r')
    return(circle_fh)

def parse_db_paths(db_list):
    name2path = {}
    for line in db_list:
        values = line.split(",")
        path = values[1].rstrip()
        name2path[values[0]] = path
    return(name2path)
                    
if __name__ == '__main__':
        '''
        This program will take a FASTA file of contigs in nucleotide form and output a table contianining circulairty and ORF information. The format of the table is as follows: contigName, length, circular, numOfORFS, numOfORFSMatchingViralFamily, bestViralFamilyClassification, numOFIntegrase, numOfACLAMEPhageParts, numOfVFDB.
        '''

        args = parse_arguments()
        
        basename = os.path.splitext(args.contigFile.name)[0]        
        circle_fh = run_circular(basename, args.contigFile)

        run_glimmer(args.contigFile)
        
        ORF_file = basename + ".fastp"
        ORF_file_fh = open(ORF_file, 'r')

        viralp_Blast_fh = run_blast_viraldb(basename, ORF_file, args.ref_viral)

        cdd_fh = run_CDD(basename, ORF_file, args.ref_cdd_db)

        ref_protein_DBs = parse_db_paths(args.ref_protein_db)
        protein2fh = {}
        for name, path in ref_protein_DBs.items():
            protein2fh[name] = run_blastp_against_db(basename, name, ORF_file, path)

        nucleotide2fh = {}
        ref_nucleotide_DBs = parse_db_paths(args.ref_nucleotide_db)    
        for name, path in ref_nucleotide_DBs.items():
            nucleotide2fh[name] = run_blastn_against_db(basename, name, args.contigFile.name, path)
            
        table = extract_annotations(args.contigFile, circle_fh, ORF_file_fh, viralp_Blast_fh, protein2fh, nucleotide2fh, cdd_fh)

        if (args.outputFile): 
            output = open(args.outputFile, 'w')
        else:
            output = sys.stdout
            
        for line in table:
            output.write('\t'.join(map(str, line)))
            output.write('\n')
