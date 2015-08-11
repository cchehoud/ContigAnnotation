import ConfigParser
import sys

class Configuration:
    """ reads ini file and parses it."""

    def __init__(self, configurationFile):
        config = ConfigParser.ConfigParser()
        config.read(configurationFile)

        self.ref_viral = config.get('Taxonomy', 'viral_protein')

        self.ref_cdd_db = config.get('CDD', 'cdd')
        self.rpsbproc_ini = config.get('CDD', 'rpsbproc_ini')
        self.skip_cdd = config.get('CDD', 'skip')
        self.skip_cdd = self.skip(self.skip_cdd)

        self.ref_protein_db = config.items('ProteinDB')
        self.ref_nucleotide_db = config.items('NucleotideDB')

    def skip(self, skip_str):
        if skip_str == 'TRUE':
            return True
        elif skip_str == 'FALSE':
            return False
        else:
            print 'skip in CDD section should be equal to : "TRUE" or "FALSE"'
            sys.exit(1)
