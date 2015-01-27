#!/usr/bin/python

import argparse

def extract_family_name (filein):
    if filein is None:
        return {}
    filein.seek(0)
    contigs2familyName = {}
    for line in filein:
        values = line.split("\t")
        contig_name = values[0].split(".")[0]
        family = values[12].split("\n")[0]
        if contig_name in contigs2familyName:
            if family in contigs2familyName[contig_name]:
                contigs2familyName[contig_name][family] += 1
            else:
                contigs2familyName[contig_name][family] = 1        
        else:
            contigs2familyName[contig_name] = {}
            contigs2familyName[contig_name][family] = 1
    return(best_family_name(contigs2familyName))

def best_family_name (contigs2familyName):
    contigs2bestFamilyName = {}
    for contig, map_of_names in contigs2familyName.items():
         highest_name = ""
         highest_count = 0
         for name, counts in map_of_names.items():
             if (counts > highest_count):
                 highest_name = name
                 highest_count = counts
         contigs2bestFamilyName[contig] = highest_name
    return(contigs2bestFamilyName)
            
if __name__ == '__main__':
    '''
    This program will take BLAST output file (format 7 without comments), as outputed by blast.viral.families.sh, and find the best viral family classification. 
    '''

    parser = argparse.ArgumentParser(description='Contig Name to Best Viral Family Classification')
    parser.add_argument('-i', '--input', dest='blastFile', required=True, type=file,
                         help='Required input: BLAST file in format 7 with no # lines.')

    args = parser.parse_args()

    print "Extracting information from BLAST file:"
    contigs2familyName = extract_family_name(args.blastFile)

    print "ContigName\tBestFamilyName"
    for contig, name in contigs2familyName.items():
        print contig + "\t" + str(name)
