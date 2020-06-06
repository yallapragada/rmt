from file_tools import open_file_by_mimetype


def read_fasta(fastafn):

    """
    Read a fasta file, return a dictionary keyed by sequence name containing the respective sequences.
    The file can be bzip or gzip compressed, or uncompressed.  Compression is detected by guessing MIME type
    """

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


