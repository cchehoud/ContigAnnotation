import sys
import optparse

def find_circular_by_kmer(filein,fo, kmerMin, kmerMax, min_len):
    found = 0
    label=''
    seq=''
    for line in filein:
        if line[0] == '>':
            if len(seq) > min_len:
                for kmer in range(kmerMin, kmerMax + 1):
                    if seq[0:kmer] == seq[(len(seq)-kmer):len(seq)]:
                        found+=1
                        fo.write(label+'\n'+seq+'\n')
                        break
            label=line.strip()
            seq=''
        else:
            seq+=line.strip()
    return found

if __name__ == '__main__':
    '''
    This program will take a FASTA file and extract the circular sequences according to a repeated k-mer at the beginning and end.
    '''

    usage = '%prog {options} sequences'
    p = optparse.OptionParser(usage=usage)
    p.add_option('-i', '--input', default=None,
        help='FASTA file.')
    p.add_option('-o', '--output_fp', default=None,
        help='Output filepath. Default is stdout.')
    p.add_option('-k', '--kmerMin', default=None,
        help='min K-mer length to look for at the beginning and end of each contig.')
    p.add_option('-m', '--kmerMax', default=None,
        help='max K-mer length to look for at the beginning and end of each contig.')
    p.add_option('-l', '--min_len', default=2000,
        help='Minimum length of circular contig.')
    opts, args = p.parse_args()

    if opts.output_fp is None:
        fo=sys.stdout
    else:
        fo=open(opts.output_fp, 'w')        
    if opts.input is None:
        p.error('Must provide FASTA file')

    if opts.kmerMin != None and opts.kmerMax != None:
        filein = open(opts.input)
        num_circular = find_circular_by_kmer(filein, fo, int(opts.kmerMin),
             int(opts.kmerMax), int(opts.min_len))
        print "Found " + str(num_circular) + " circular contigs."
        filein.close()
        fo.close()
    else: p.error('Must provide kmer length')

