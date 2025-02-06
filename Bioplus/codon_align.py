from Bio import SeqIO
import subprocess
import argparse
import tempfile
from io import StringIO
import os
import sys


def align_amino(records, bin='mafft'):
    handle = tempfile.NamedTemporaryFile(delete=False)
    for record in records:
        aaseq = record.seq.translate()
        line = f">{record.name}\n{aaseq}\n"
        handle.write(line.encode('utf-8'))
    handle.close()

    stdout = subprocess.check_output([bin, '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    os.remove(handle.name)  # delete temp file
    return SeqIO.parse(StringIO(stdout), "fasta")


def apply_aln(records, aln):
    """
    Apply the alignment of translated amino acid sequences to the nucleotide
    sequences to yield a codon-aware alignment.
    :param records:  list, SeqRecord objects representing unaligned nucleotide
                     sequences
    :param aln:  SeqIO iterable object from align_amino()
    :yield:  tuple, sequence name and codon-aligned nucleotide sequence
    """
    for i, aaseq in enumerate(aln):
        nucseq = records[i].seq  # original nucleotide sequence
        newseq = ''
        pos = 0
        for aa in aaseq:
            if aa == '-':
                newseq += '---'
                continue
            newseq += nucseq[pos:(pos+3)]
            pos += 3
        yield records[i].name, newseq


if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser("Codon-aware alignment")
    parser.add_argument(
        "infile", type=argparse.FileType('r'),
        help="input, path to file containing homologous codon sequences, i.e., "
             "nucleotide sequences that are copies of the same open reading frame."
    )
    parser.add_argument(
        "-f", "--format", type=str, default='fasta',
        help="option, format of input file (must be supported by Bio.SeqIO)."
    )
    parser.add_argument(
        "-b", "-bin", type=str, default="mafft",
        help="option, path to MAFFT binary file."
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="option, path to write FASTA output. Default is stdout."
    )
    args = parser.parse_args()

    # load records into memory because we need to iterate twice
    records = list(SeqIO.parse(args.infile, args.format))
    aln = align_amino(records)
    for label, codseq in apply_aln(records, aln):
        args.outfile.write(f">{label}\n{codseq}\n")
