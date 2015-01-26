import unittest

from .. import contigs2length

class Test_fasta_length(unittest.TestCase):
    """ should be able to get length of read/contig for different versions
        of fasta format: seq on one line, seq on separate lines
    """

    def setUp(self):
        self.one_line_seq = open("tests/one_line_seq.fa")
        self.repeated_id = open("tests/one_line_seq_id_repeated.fa")
        self.multiline_seq = open("tests/multiline_seq.fa")

    def test_number_of_reads(self):
        self.assertEquals(
            len(contigs2length.extract_name_length(self.one_line_seq)), 3)
        self.assertEquals(
            len(contigs2length.extract_name_length(self.multiline_seq)), 2)

    def test_oneline_fasta_counts(self):
        name_length = contigs2length.extract_name_length(self.one_line_seq)
        self.assertEquals(name_length["length_116"], 116)
        self.assertEquals(name_length["length_15"], 15)
        self.assertEquals(name_length["length_200"], 200)

    def test_multiline_fasta_counts(self):
        name_length = contigs2length.extract_name_length(self.multiline_seq)
        self.assertEquals(name_length["length_116"], 116)
        self.assertEquals(name_length["length_200"], 200)

    def test_all_ids_are_unique(self):
        self.assertRaises(ValueError, contigs2length.extract_name_length, self.repeated_id)
