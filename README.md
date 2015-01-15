# Contig Annotation

For each nucleotide sequence collects annotation by searching nucleotide, protein and conserved domain
databases. Each nucleotide sequence translated to protein sequence after ORFs are predicted.


Here is example command to run the script:

```bash
python contigAnnotation.py -i ex_contig.fa -o ex_contig_annotation.dat -p proteinRefDB.txt\ 
-n nucleotideRefDB.txt -r /media/THING1/sminot/timecourse/4AnnotateContigs/4.12TaxonomicFamily/4.12.1ViralFamilyProteinsDB/\
-c /media/THING1/dryga/PhageDynamics/CDD/cdd/little_endian 
```

The command-line options are not stable and could be changed.

## Input 

Input files for the script are:
* fasta (nucleotide) file with contigs
* file with paths to blastp protein databases
* blastp protein databases 
* file with paths to blastn nucleotide databases
* blastn nucleotide databases
* file with protein databases for viral families

Format for files with blast DBs, both (protein and nucleotide) are
tab-delimited with 2 columns: name of the database and path to it.
See example files in proteinRefDB.txt and nucleotideRefDB.txt.

CDD utility needs file `rpsbproc.ini`, which should be located in the current directory and point to
correct location, where CDD files can be found.

## Output

Creates a tab-delimited file, where first line is a header line
and rest lines have annotation for each sequence from FASTA input file. 

Provides information about length, ORFs, viral DB matches, number of hits to each protein DB, top hit for each nucleotide DB, and number of hits for CDD DB.
 
The current format is:
```
contig_name contig_length, number of ORFs,... 
```

## Dependencies

Script uses several programs and databases to add annotation to sequences.

### Programs

* blastn
* blastp
* rpsblast
* rpsbproc
* glimmer

### Databases

Blast databases are required for annotating sequences with blast hits.
CDD database and additional files reqauired by `rpsbproc.ini`

