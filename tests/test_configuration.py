
import unittest
import StringIO
from .. import configuration

class Test_configuration(unittest.TestCase):
    
    def setUp(self):
        valid_config_file = "tests/valid_config_file.ini"
        self.configuration = configuration.Configuration(valid_config_file)

    def test_ref_viral(self):
        ref_viral =  '/media/THING1/sminot/timecourse/4AnnotateContigs/4.12TaxonomicFamily/4.12.1ViralFamilyProteinsDB/'
        self.assertEquals(self.configuration.ref_viral, ref_viral) 

    def test_ref_cdd(self):
        cdd = '/media/THING1/dryga/PhageDynamics/CDD/cdd/little_endian' 
        rpsbproc_ini = './rpsbproc.ini'
        self.assertEquals(self.configuration.ref_cdd_db, cdd) 
        self.assertEquals(self.configuration.rpsbproc_ini, rpsbproc_ini) 

    def test_protein_db(self):
        self.assertEquals(len(self.configuration.ref_protein_db), 3)
        self.assertEquals(self.configuration.ref_protein_db[1], 
            ("aclame",  '/media/THING1/local/genomeIndexes/blast/ACLAME/aclame_proteins_viruses_prophages_0.4.fasta'))

        self.assertEquals(self.configuration.ref_protein_db[2], 
            ("vfdb",  '/media/THING1/local/genomeIndexes/blast/VFDB/VFs.faa'))
    
    def test_nucleotide_db(self):
        self.assertEquals(len(self.configuration.ref_nucleotide_db), 3)
        self.assertEquals(self.configuration.ref_nucleotide_db[0], 
             ('viral', '/home/rohinis/viral_blastdb/viral.1.1.genomic.fna'))
        self.assertEquals(self.configuration.ref_nucleotide_db[1], 
             ('nt', '/media/THING1/local/genomeIndexes/blast_nt/nt'))
        self.assertEquals(self.configuration.ref_nucleotide_db[2], 
             ('bacteria', '/media/THING1/local/genomeIndexes/blast/BacterialGenomes/ncbi_bacteria.fa'))
        

        
        
