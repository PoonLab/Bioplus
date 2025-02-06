import math
from Bio import AlignIO
import argparse
import csv
import sys

description = """\
Calculate a k-mer distance (Euclidean, spectrum kernel), p-distance or 
Jukes-Cantor corrected distance for an alignment of nucleotide sequences."
"""


def kmer(seq, k):
    """
    Calculate word counts for words of length k in sequence seq.
    :param seq: str, nucleotide sequence
    :param k: int, word length
    :return: dict, counts keyed by word - absent entries imply zero
    """
    d = {}
    for i in range(k, len(seq)):
        word = seq[(i-k):i]
        if word not in d:
            d.update({word: 0})
        d[word] += 1
    return d


def kdist(s1, s2, cache, kernel=False):
    """
    Calculate Euclidean or spectrum kernel distance.
    Leslie C, Eskin E, Noble WS. The spectrum kernel: a string kernel for SVM
    protein classification. Pac Symp Biocomput. 2002:564-75. PMID: 11928508.

    :param s1:  str, sequence 1
    :param s2:  str, sequence 2
    :param cache:  dict, precomputed word counts
    :param kernel:  bool, if True return spectrum kernel
    :return:  float, pairwise distance
    """
    k1 = cache[s1]['k']
    k2 = cache[s2]['k']
    res = 0.
    words = set(k1.keys()).union(set(k2.keys()))
    for word in words:
        c1 = k1.get(word, 0)
        c2 = k2.get(word, 0)
        if kernel:
            res += c1 * c2  # dot product
        else:
            res += (c1 - c2)**2  # Euclidean

    if kernel:
        return 1 - res / math.sqrt(cache[s1]['norm'] * cache[s2]['norm'])
    return math.sqrt(res)


def pdist(s1, s2, corrected=False):
    """ Calculate the p-distance or Jukes-Cantor (corrected) """
    ndiff = 0
    for i, nt1 in enumerate(s1):
        nt2 = s2[i]
        if nt1 != nt2 and nt1 not in '-N' and nt2 not in '-N':
            ndiff += 1
    p = ndiff / len(s1)
    if corrected:
        return abs(-3 / 4 * math.log(1 - 4 / 3 * p))
    return p


if __name__ == "__main__":
    # command-line interace
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "infile", type=argparse.FileType('r'),
        help="input, FASTA file with aligned nuc sequences"
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="option, path to write CSV output. Default is stdout."
    )
    parser.add_argument(
        "-d", "--dist", choices=["euc", "spec", "pdist", "jc"],
        help="Specify which distance to calculate: Euclidean (euc), spectrum "
             "kernel (spec), p-distance (pdist), or Jukes-Cantor (jc, default)."
    )
    parser.add_argument(
        "-k", type=int, default=5,
        help="option, k-mer length (default 5 nt, applies to euc/spec)"
    )
    parser.add_argument(
        "-f", "--format", default="fasta",
        help="format of sequence file format, passed to SeqIO.parse"
    )
    parser.add_argument(
        "-n", action="store_true",
        help="option, write numerical indices instead of sequence headers. "
             "Index reflects ordering of records in the input file."
    )
    args = parser.parse_args()

    # extract sequence headers (labels)
    records = AlignIO.read(args.infile, args.format)
    seqs, labels = [], []
    for record in records:
        seqs.append(str(record.seq))
        labels.append(record.name)

    # if computing a kmer distance, extract words first
    cache = {}
    if args.dist in ["euc", "spec"]:
        for s in seqs:
            if s not in cache:  # skip duplicate sequences
                kdict = kmer(s, args.k)
                norm = sum([count*count for w, count in kdict.items()])
                cache.update({s: {"k": kdict, "norm": norm}})

    # calculate distances and write outputs
    writer = csv.writer(args.outfile)
    writer.writerow(['seq1', 'seq2', 'dist'])
    for i in range(0, len(seqs)-1):
        for j in range(i+1, len(seqs)):
            # calculate pairwise distance
            if args.dist in ["euc", "spec"]:
                d = kdist(seqs[i], seqs[j], cache=cache,
                              kernel=(args.dist=="spec"))
            else:
                d = pdist(seqs[i], seqs[j], corrected=(args.dist=='jc'))

            writer.writerow([i, j, d] if args.n else [labels[i], labels[j], d])
