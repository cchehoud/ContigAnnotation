#!/usr/bin/python

import argparse

def extract_domain_counts(filein, contig_cdd_accession_fh):
    """ extract counts and save contig-cdd_id table into file."""
    if filein is None:
        return {}
    filein.seek(0)
    contigs_2_domain_count = {}
    for line in filein:
        if (line.startswith("QUERY")):
            contig = extract_contig_name(line)
            count = 0
            line = next(filein)
            if (line.startswith("DOMAINS")):
                line = next(filein)
                while not (line.startswith("ENDDOMAINS")):
                    cdd_accession = parse_accession(line)
                    contig_cdd_accession_fh.write(
                        contig + "\t" + cdd_accession + "\n")
                    count+=1
                    line = next(filein)
            if contig in contigs_2_domain_count:
                contigs_2_domain_count[contig]+= count
            else:
                contigs_2_domain_count[contig] = count
    return(contigs_2_domain_count)

def parse_accession(line):
    accesion_index = 8
    row = line.split()
    if len(row) < accesion_index + 1:
        raise ValueError("cannot process rpsbproc output file.")
    return row[accesion_index]

def extract_contig_name(line):
    values = line.split("\t")
    return(values[4].split(".")[0])

if __name__ == '__main__':
    '''
    This program will take a rpsbproc output file and return a map of contig to number of domains.
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Number of domains.')
    parser.add_argument('-i', '--input', dest='rpsbprocFile', required=True, type=file,
                         help='Required input: rpsbproc file.')

    args = parser.parse_args()

    print "Extracting information from rpsbproc file:"
    contigs2domaincount = extract_domain_counts(args.rpsbprocFile, open("contig_cdd_accession.tsv", 'w'))

    print "ContigName\tNumberofDomains"
    for contig, nDomains in contigs2domaincount.items():
        print contig + "\t" + str(nDomains)
