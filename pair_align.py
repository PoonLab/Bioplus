from io import StringIO
from Bio import SeqIO
import argparse
import sys
import tempfile
import subprocess

description = """\
Use MAFFT to align sequences to a reference sequence, discarding any 
insertions.  This can also be used to trim sequences to a specific 
region of the genome, such as a reference gene.
"""


def mafft(query, ref, binpath="mafft"):
    """
    :param query:  str, sequence to align to reference
    :param ref:  str, reference sequence
    :param binpath:  str, path to MAFFT executable file
    :return:  dict
    """
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = f'>ref\n{ref}\n>query\n{query}\n'
    handle.write(s.encode('utf-8'))
    handle.close()

    # call MAFFT on temporary file
    stdout = subprocess.check_output([binpath, '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    result = list(SeqIO.parse(StringIO(stdout), "fasta"))
    aref = str(result[0].seq)
    aquery = str(result[1].seq)

    newseq = ""
    for i, aa in enumerate(aref):
        if aa == '-':
            continue
        newseq += aquery[i]
    return newseq


if __name__ == "__main__":
    # command-line interface
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "infile", type=argparse.FileType('r'),
        help="input, file with unaligned sequences"
    )
    parser.add_argument(
        "ref", type=argparse.FileType('r'),
        help="input, file containing reference to align sequences to."
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="option, path to write CSV output. Default is stdout."
    )
    parser.add_argument(
        "-b", "--binpath", default="mafft",
        help="Path to MAFFT executable file.  Defaults to 'mafft'."
    )
    parser.add_argument(
        "-f", "--format", default="fasta",
        help="format of sequence file format, passed to SeqIO.parse"
    )
    args = parser.parse_args()

    # read reference sequence
    reference = SeqIO.read(args.ref, args.format)
    ref = str(reference.seq)

    # extract sequence headers (labels)
    records = SeqIO.parse(args.infile, args.format)
    for record in records:
        sys.stderr.write(record.name+'\n')
        sys.stderr.flush()

        query = str(record.seq)
        aligned = mafft(query=query, ref=ref, binpath=args.binpath)
        args.outfile.write(f">{record.description}\n{aligned}\n")
