import os
import unittest
import tempfile

from .. import fasta_filter
from .. import contigs2length

class Test_Fasta_Filtering(unittest.TestCase):
    
    def setUp(self):
        self.fasta = "tests/one_line_seq.fa"
        self.output = "tests/filtered_fasta.fa"

    def tearDown(self):
        if os.path.isfile(self.output):
            os.remove(self.output)

    def test_raises_if_no_file(self):
        self.assertRaises(IOError, fasta_filter.seq_length_greater, "NonExistingFile.fasta", self.output, 42)
        self.assertTrue(not os.path.isfile(self.output))

    def test_all_reads_survive(self):
        self.assertEquals(fasta_filter.seq_length_greater(self.fasta, self.output, 0), 3)
        self.assertEquals(len(contigs2length.extract_name_length(open(self.output))), 3) # from file

    def test_some_reads_survive(self):
        self.assertEquals(fasta_filter.seq_length_greater(self.fasta, self.output, 100), 2)
        self.assertEquals(len(contigs2length.extract_name_length(open(self.output))), 2) # from file

    def test_one_read_survive(self):
        self.assertEquals(fasta_filter.seq_length_greater(self.fasta, self.output, 200), 1)
        self.assertEquals(len(contigs2length.extract_name_length(open(self.output))), 1) # from file

    def test_none_of_reads_survive(self):
        self.assertEquals(fasta_filter.seq_length_greater(self.fasta, self.output, 500), 0)
        self.assertTrue(not os.path.isfile(self.output)) # nothing written to file because all read are smaller

