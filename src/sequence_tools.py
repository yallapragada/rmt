from file_tools import open_file_by_mimetype

def rev_comp(seq, nowarn = 0):
    '''Accepts a sequence, forces them to capitals, and returns the reverse complement; accepted 
    characters in the sequence are A, T, C, G, N
    '''
    rcseq = seq.upper()

    rcseq = rcseq.replace('A', 'z')
    rcseq = rcseq.replace('T', 'A')
    rcseq = rcseq.replace('z', 'T')

    rcseq = rcseq.replace('C', 'z')
    rcseq = rcseq.replace('G', 'C')
    rcseq = rcseq.replace('z', 'G')

    return rcseq[::-1]


def read_fasta(fastafn):
    '''Read a fasta file, return a dictionary keyed by sequence name containing the respective sequences.

    The file can be bzip or gzip compressed, or uncompressed.  Compression is detected by guessing MIME type.'''

    fastafh = open_file_by_mimetype(fastafn, 'r')

    # parse the FASTA file
    seqs = {}
    rows = []
    seqname = ''
    for line in fastafh:
        if line.startswith('>'):
            if seqname != '':
                seqs[seqname] = ''.join(rows)
                seqname = line.strip()[1:]
                rows = []
            else:
                seqname = line.strip()[1:]
                rows = []
        else:
            rows.append(line.strip())
    seqs[seqname] = ''.join(rows)

    return seqs


