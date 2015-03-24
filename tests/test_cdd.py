import os
import unittest
import tempfile
import csv

from .. import contigs2CDD

def check_pair_exists(pair, list_of_pairs):
    return pair in list_of_pairs

class Test_Fasta_Filtering(unittest.TestCase):
    
    def setUp(self):
        self.cdd_input = open("tests/cdd.txt")
        self.cdd_output = open("tests/cdd_accession.tsv", 'w+b')

    def tearDown(self):
        if os.path.isfile(self.cdd_output.name):
            os.remove(self.cdd_output.name)
        self.cdd_input.close()

    def test_count_correct(self):
        self.assertEquals(len(contigs2CDD.extract_domain_counts(
            self.cdd_input, self.cdd_output)), 3)

    def test_cdd_accession_file_created(self):
        contigs2CDD.extract_domain_counts(self.cdd_input, self.cdd_output)
        self.assertTrue(os.path.exists(self.cdd_output.name))

    def test_cdd_accession_file_with_correct_number_of_pairs(self):
        contigs2CDD.extract_domain_counts(self.cdd_input, self.cdd_output)
        n_lines = 0 
        self.cdd_output.seek(0)
        for line in self.cdd_output:
            n_lines += 1
        self.assertEquals(n_lines, 7) 

    def test_cdd_accession_file_number_of_rows(self):
        contigs2CDD.extract_domain_counts(self.cdd_input, self.cdd_output)
        self.cdd_output.seek(0)
        reader = csv.reader(self.cdd_output, delimiter="\t")
        for row in reader:
            self.assertEquals(len(row), 2) 

    def test_cdd_accession_pairs_are_correct(self):
        contigs2CDD.extract_domain_counts(self.cdd_input, self.cdd_output)
        self.cdd_output.seek(0)
        reader = csv.reader(self.cdd_output, delimiter="\t")
        contig_accession = {}
        for row in reader:
            contig_accession[row[0]] = row[1]
        self.assertTrue("contig-95_2" in contig_accession)
        self.assertTrue("contig-95_4" in contig_accession)
        self.assertTrue("contig-NO-HITS" not in contig_accession)

    def test_cdd_accession_pairs_number(self):
        contigs2CDD.extract_domain_counts(self.cdd_input, self.cdd_output)
        self.cdd_output.seek(0)
        reader = csv.reader(self.cdd_output, delimiter="\t")
        contig_accession = {}
        for row in reader:
            contig = row[0]
            if contig not in contig_accession:
                contig_accession[contig] = 1
            else:
                contig_accession[contig] += 1
        self.assertEquals(contig_accession["contig-95_2"], 3)
        self.assertEquals(contig_accession["contig-95_4"], 4)

    def test_cdd_accession_pairs_values(self):
        contigs2CDD.extract_domain_counts(self.cdd_input, self.cdd_output)
        self.cdd_output.seek(0)
        reader = csv.reader(self.cdd_output, delimiter="\t")
        contig_accession = []
        for row in reader:
            contig_accession.append(tuple(row))
        self.assertTrue(check_pair_exists( ("contig-95_2", "pfam02956"), contig_accession))
        self.assertTrue(check_pair_exists( ("contig-95_2", "pfam10145"), contig_accession))
        self.assertTrue(check_pair_exists( ("contig-95_2", "COG5280"  ), contig_accession))

        self.assertFalse(check_pair_exists( ("contig-95_2", "COG3500" ), contig_accession))

        self.assertFalse(check_pair_exists( ("contig-95_4", "pfam02956"), contig_accession))
        self.assertFalse(check_pair_exists( ("contig-95_4", "COG5280"  ), contig_accession))

        self.assertTrue(check_pair_exists( ("contig-95_4", "COG3500"  ), contig_accession))
        self.assertTrue(check_pair_exists( ("contig-95_4", "pfam10145"), contig_accession))
        self.assertTrue(check_pair_exists( ("contig-95_4", "cl00465"  ), contig_accession))
        self.assertTrue(check_pair_exists( ("contig-95_4", "cl19582"  ), contig_accession))

