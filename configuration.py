import ConfigParser

class Configuration:
    """ reads ini file and parses it 
        
    """

    def __init__(self, configurationFile):
        config = ConfigParser.ConfigParser()
        config.read(configurationFile)

        self.ref_viral = config.get('Taxonomy', 'viral_protein')

        self.ref_cdd_db = config.get('CDD', 'cdd')
        self.rpsbproc_ini = config.get('CDD', 'rpsbproc_ini')

        self.ref_protein_db = config.items('ProteinDB')
        self.ref_nucleotide_db = config.items('NucleotideDB')
