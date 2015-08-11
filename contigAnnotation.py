#!/usr/bin/python

import sys
import argparse
import subprocess
import os
import shutil

import circular
import contigs2length
import contigs2circular
import contigs2ORFs
import contig2ViralFamily
import contigs2count
import contigs2topMatch
import contigs2CDD
import configuration 
import fasta_filter
import viral_host

def run_CDD(basename, ORF_file, ref_CDD_DB, rpsbproc_ini):
    CDD_xml_file = basename + "_cdd.xml"
    print "Running CDD search:"
    command = "rpsblast -query " + ORF_file + " -evalue 0.01 -seg no -outfmt 5 -db " + ref_CDD_DB + " -num_threads 4 -out " + CDD_xml_file
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError:
        print "CDD annotation(rpsblast) does not return any hits."
        return None
    return run_rpsbproc(CDD_xml_file, rpsbproc_ini)

def run_rpsbproc(CDD_xml_file, rpsbproc_ini):
    """ run utility for blast output postprocesing.
       
        Temporary workaround for cmd-line processing of rpsbproc:
        input file and output file are in the current directory. 
    """
    current_path = os.getcwd()
    try:
        (path_to_cdd, filename_extension) = os.path.split(CDD_xml_file)
        CDD_xml_file = filename_extension
        CDD_out_file = os.path.splitext(filename_extension)[0] + ".txt" 
        (rpsbproc_ini_path, rpsbproc_ini_filename) =  os.path.split(rpsbproc_ini)
        command2 = "rpsbproc -i " + CDD_xml_file + " -o " + CDD_out_file + " -c " + rpsbproc_ini_filename
        try:
            rpsbproc_copy_to = os.path.join(path_to_cdd, rpsbproc_ini_filename)
            shutil.copyfile(rpsbproc_ini, rpsbproc_copy_to)
        except shutil.Error:
            # if cannot copy file it is in the same folder already
            pass
        if path_to_cdd != '': # empty when current file is in current dir and only name of the file is given
            os.chdir(path_to_cdd) 
        subprocess.check_call(command2, shell=True)
        return open(CDD_out_file, 'r')
    except subprocess.CalledProcessError:
        print "CDD annotation(rpsbproc) does not extract any CDD domains."
        return None
    finally:
        os.chdir(current_path)
        
        
def run_glimmer(contig_file):
    """ run ORF finder and creates file with predicted protein seq.

        Return: True if success
                False if cannot run glimmer of get protein seqs.
    """
    print "Running Glimmer:"
    try:
        subprocess.check_call(["bash", "glimmer-wrapper.sh", contig_file])
        return True
    except subprocess.CalledProcessError:
        return False
 
def run_blast_viraldb(basename, ORF_file, ref_viralDB):
    viralp_Blast_file = basename + "_viralp_blastout.txt"
    print "Running BLAST on viral proteins db:"
    try:
        subprocess.check_call(["bash", "blast.viral.families.sh", ORF_file, viralp_Blast_file, ref_viralDB])
        return open(viralp_Blast_file, 'r')
    except subprocess.CalledProcessError:
        print "blast against viral families does not return any hits."
        return None
    
def run_blastp_against_db(basename, name, ORF_file, path):
    output_Blast_file = basename + "_" + name + "_blastout.txt"
    print "Running BLAST on " + name + " db:"
    command = "blastp -query " + ORF_file + " -db " + path + " -outfmt '6 qseqid' -evalue 1e-10 -num_threads 4 -max_target_seqs 1 -out " + output_Blast_file
    try:
        subprocess.check_call(command, shell=True)
        return open(output_Blast_file, 'r')
    except subprocess.CalledProcessError:
        print "Protein BLAST(blastp) failed for DB: " + path
        return None
    
def run_blastn_against_db(basename, name, contig_file, path):
    output_Blast_file = basename + "_" + name + "_blastout.txt"
    print "Running BLAST on " + name + " db:"
    command = "blastn -query " + contig_file + " -db " + path + " -outfmt '6 qseqid stitle sseqid length pident' -num_threads 4 -evalue 1e-10 -max_target_seqs 1 -out " + output_Blast_file
    try:
        subprocess.check_call(command, shell=True)
        return open(output_Blast_file, 'r')
    except subprocess.CalledProcessError:
        print "Nucleotide BLAST(blastn) failed for DB: " + path
        return None
         
def extract_annotations(basename, contig_file, circle_file_fh, orf_file_fh, viralp_Blast_fh, protein2fh, nucleotide2fh, cdd_fh): 
    contig_file_fh = open(contig_file)
    contig_cdd_accession = basename + "_cdd_accession.tsv"
    print "contig cdd_id table will be written to : " + contig_cdd_accession
    c2length = contigs2length.extract_name_length(contig_file_fh)
    c_circular = contigs2circular.extract_circularity(circle_file_fh)
    c2ORFs = contigs2ORFs.extract_ORF_counts(orf_file_fh)
    c2viralORFs = contigs2count.extract_counts(viralp_Blast_fh)
    c2familyName = contig2ViralFamily.extract_family_name(viralp_Blast_fh)
    c2domains = contigs2CDD.extract_domain_counts(cdd_fh, open(contig_cdd_accession, 'w'))
    
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
    header += proteinDBmatches.keys() + nucleotideDBmatches.keys() + ["host"]
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
        host = viral_host.extract_putative_host(nucleotideDBmatches.keys(), nuc_list)
        row = [contig, length, circular, nORFs, nViralORFs, Family, num_domains]
        row += protien_list + nuc_list + [host]
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
    parser.add_argument('-d' '--databases', dest='databases', required=True, type=str,
                        help='Required input: configuration file(INI) with reference names and paths.')
    parser.add_argument('-l' '--min_length', dest='contig_min_length', required=False, type=int,
                        help='only contigs >= given length are annotated.')
    return(parser.parse_args())

def run_circular(basename, contig_file):
    print "Running Circularity check:"
    circle_file = basename + "_circularity.fasta"
    circle_fh = open(circle_file, 'w')
    kMerMin = 10
    kMerMax = 1000
    min_len = 3500
    circular.find_circular_by_kmer(open(contig_file), circle_fh, kMerMin, kMerMax, min_len)
    circle_fh.close()
    circle_fh = open(circle_file, 'r')
    return(circle_fh)

def run_protein_searches(basename, configuration):
        ORF_file = basename + ".fastp"
        ORF_file_fh = open(ORF_file, 'r')

        viralp_Blast_fh = run_blast_viraldb(basename, ORF_file, configuration.ref_viral)

        cdd_fh = None
        if configuration.skip_cdd:
            print "CDD searches will NOT be done."
        else:
            print "CDD searches will be done."
            cdd_fh = run_CDD(basename, ORF_file, configuration.ref_cdd_db, configuration.rpsbproc_ini)

        protein2fh = {}
        for (name, path) in configuration.ref_protein_db:
            protein2fh[name] = run_blastp_against_db(basename, name, ORF_file, path)
        return (ORF_file_fh, viralp_Blast_fh, protein2fh, cdd_fh)

def filter_contigs(contig_file, min_length):
    """ created file and names for filtered contig

        Exit if none of the contigs survive filtering.
    
        Return:
            (basename, contig_name)
    """
    basename_original = os.path.splitext(contig_file)[0]        
    basename_filtered = basename_original + "_min_len_" + str(min_length)
    contigs_filtered = basename_filtered + ".fasta"
    num_reads = fasta_filter.seq_length_greater(contig_file, contigs_filtered, min_length)
    if num_reads == 0:
        print "All of the contigs are smaller than " + str(min_length) + " nucleotide."
        print "NOTHING is going to be annotated. (Decrease length or drop the argument for annotation.)"
        sys.exit(0)
    return (basename_filtered, contigs_filtered)
    
                    
if __name__ == '__main__':
        '''
        This program will take a FASTA file of contigs in nucleotide form and output a table contianining circulairty and ORF information. The format of the table is as follows: contigName, length, circular, numOfORFS, numOfORFSMatchingViralFamily, bestViralFamilyClassification, numOFIntegrase, numOfACLAMEPhageParts, numOfVFDB.
        '''

        args = parse_arguments()
        configuration = configuration.Configuration(args.databases)
        
        if args.contig_min_length:
            (basename, contigs) = filter_contigs(args.contigFile.name, args.contig_min_length)
        else:
            contigs = args.contigFile.name
            basename = os.path.splitext(contigs)[0]        

        circle_fh = run_circular(basename, contigs)

        is_protein_seq_created = run_glimmer(contigs)
        ORF_file_fh = None
        viralp_Blast_fh = None
        protein2fh = {}
        cdd_fh = None
        if not is_protein_seq_created:
            print "Glimmer FAILED: all searches based on protein seq will be SKIPPED."
        else:
            print "Glimmer succeed: all searches based on protein seq will be RUN."
            (ORF_file_fh, viralp_Blast_fh, protein2fh, cdd_fh) =  run_protein_searches(basename, configuration)

        nucleotide2fh = {}
        ref_nucleotide_DBs = configuration.ref_nucleotide_db    
        for name, path in ref_nucleotide_DBs:
            nucleotide2fh[name] = run_blastn_against_db(basename, name, contigs, path)
            
        table = extract_annotations(
            basename, contigs, circle_fh, ORF_file_fh, viralp_Blast_fh, protein2fh, nucleotide2fh, cdd_fh)

        if (args.outputFile): 
            output = open(args.outputFile, 'w')
        else:
            output = sys.stdout
            
        for line in table:
            output.write('\t'.join(map(str, line)))
            output.write('\n')

