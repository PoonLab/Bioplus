import random
import argparse


description = """
Random permutation of FASTA file.
This script reads a FASTA-formatted file containing aligned sequences and
applies a random permutation of alignment columns (nucleotide sites),
writing the result to a new FASTA file.  A random permutation changes the
order of a sequence.  For example, (1,2,3,4,5) may become (4,2,1,3,5).
This has the benefit of providing some data security (by scrambling the
original sequences) while preserving information about evolutionary
relationships.
"""


def parse_fasta(handle):
    res = []
    header = None
    sequence = ''
    for line in handle:
        if line.startswith(">"):
            if len(sequence) > 0:
                res.append([header, sequence])
                sequence = ''   # reset containers
            header = line.strip('>\n')
        else:
            sequence += line.strip()
    res.append([header, sequence])
    return res


def transpose_fasta(fasta):
    """ Convert list of sequences to list of columns """
    n_columns = len(fasta[0][1])
    columns = ['' for i in range(n_columns)]
    for _, seq in fasta:
        for i, nt in enumerate(seq):
            columns[i] += nt
    return columns


def untranspose_fasta(columns):
    """ Convert columns back to a list of sequences """
    nseq = len(columns[0])
    seqs = ['' for i in range(nseq)]
    for col in columns:
        for i, nt in enumerate(col):
            seqs[i] += nt
    return seqs


if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("input", type=argparse.FileType('r'),
                        help="input, FASTA-formatted file of aligned sequences")
    parser.add_argument("output", type=argparse.FileType('w'),
                        help="output, FASTA-formatted file of permuted sequences")
    args = parser.parse_args()

    # import sequences from FASTA file
    fasta = parse_fasta(args.input)
    headers = [h for h, s in fasta]  # save headers in a list

    # do random permutation of columns in sequence alignment
    columns = transpose_fasta(fasta)
    random.shuffle(columns)
    pseqs = untranspose_fasta(columns)
    pfasta = zip(headers, pseqs)

    # write result to another file
    for header, seq in pfasta:
        args.output.write(f">{header}\n{seq}\n")
    args.output.close()
