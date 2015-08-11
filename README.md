# Contig Annotation

For each nucleotide sequence collects annotation by searching nucleotide, protein and conserved domaini(CDD)
databases. Each nucleotide sequence translated to protein sequence after ORFs are predicted.

Here is example command to run the script:

```bash
python contigAnnotation.py -i ex_contig.fa -o ex_contig_annotation.dat -d databases.ini
```

## Input 

Input files for the script are:

* fasta (nucleotide) file with contigs
* configuration file with paths to databases

Configuration file in INI Format:

```ini
[Taxonomy]
viral_protein = /media/THING1/sminot/timecourse/4AnnotateContigs/4.12TaxonomicFamily/4.12.1ViralFamilyProteinsDB/

[CDD]
skip=TRUE
cdd = /media/THING1/dryga/PhageDynamics/CDD/cdd/little_endian 
rpsbproc_ini = ./rpsbproc.ini

[ProteinDB]
integrase = /media/THING1/local/genomeIndexes/blast/UniprotPhageIntegrase.fasta
aclame = /media/THING1/local/genomeIndexes/blast/ACLAME/aclame_proteins_viruses_prophages_0.4.fasta
vfdb = /media/THING1/local/genomeIndexes/blast/VFDB/VFs.faa

[NucleotideDB]
viral = /home/rohinis/viral_blastdb/viral.1.1.genomic.fna
nt = /media/THING1/local/genomeIndexes/blast_nt/nt
bacteria = /media/THING1/local/genomeIndexes/blast/BacterialGenomes/ncbi_bacteria.fa
```

Configuration file has 4 sections: Taxonomy, CDD, ProteinDB, NucleotideDB.

Taxonomy section has only one key/value pair, key should be `viral\_protein`
and value is path to blastp protein database for viral family.

CDD section has 3 key/value pairs, 1st key is `cdd` and value is path to CDD database,
2nd key is `rpsbproc\_ini` and value is path to init file required for rpsbproc utility. 
CDD searches can be disabled to save running time:
```
skip=TRUE 
```
to actually run CDD use:
```
skip=FALSE 
```


ProteinDB section has arbitrary number of key/value pairs, each key is name of
protein blast db and value is path to the DB.

NucleotideDB section has arbitrary number of key/value pairs, each key is name of
nucleotide blast db and value is path to the DB.

CDD utility needs file `rpsbproc.ini`.

## Output

Creates a tab-delimited file, where first line is a header line
and rest of the lines have annotation for each sequence from FASTA input file. 

Provides information about length, cirlularity, number of predicted ORFs, number of
ORFs that matches viral protein db, putative viral family based on viral matches, number of CDD domains,
number of hits to each protein DB(given in ProteinDB section of INI file), 
and top hit for each nucleotide DB(given in NucleotideDB section of INI file).
 
The current format is:
```
contig_name	length	circular	nORFs	nViralORFs	family	nDomains	integrase   ... 	aclame	viral	...	bacteria
```
Where '...' represents all the databases that are in ProteinDB(or NucleotideDB) section of INI file.

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
CDD database and additional files required by `rpsbproc.ini`.

