import os.path

from Bio import SeqIO

def seq_length_greater(input_fasta, output_fasta, minimum_length):
    """ keep all reads/contigs from input file that have at least( >= ) minimum_length nucleotides. 

        If none of the reads survive the output file is not written,
        so if return value is 0 output_fasta will not be created.

        Raises:
            IOError if input_file does not exist
        Return:
            number of surviving reads/contigs
        Side effect:
            surviving reads are written to output_fasta 
    """
    if (not os.path.isfile(input_fasta)):
        raise IOError("FASTA input file does not exist")
    number_reads = 0 # number of surviving reads
    output_handle = open(output_fasta, "w")
    for record in SeqIO.parse(input_fasta, "fasta"):
        if len(record) >= minimum_length:
            number_reads += 1
            SeqIO.write(record, output_handle, "fasta")        
    output_handle.close()
    if number_reads == 0:
        os.remove(output_fasta) # no use of empty file
    return number_reads
