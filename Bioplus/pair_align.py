from io import StringIO
from Bio import SeqIO
import argparse
import sys
import tempfile
import subprocess

try:
    from mpi4py import MPI
    nprocs = MPI.COMM_WORLD.Get_size()
    my_rank = MPI.COMM_WORLD.Get_rank()
except ModuleNotFoundError:
    sys.stderr.write("Running in serial mode...\n")
    nprocs = 1
    my_rank = 0

description = """\
Use MAFFT to align sequences to a reference sequence, discarding any 
insertions.  This can also be used to trim sequences to a specific 
region of the genome, such as a reference gene.

To run in parallel, you must have mpi4py installed and call this script
using mpirun, i.e.: 
  mpirun -np <number of processes> python3 pair_align [...]
The outputs will be written to separate files with integer suffixes.
You can concatenate these outputs into a single FASTA file with a `cat` 
command.
"""


def mafft(query, ref, binpath="mafft"):
    """
    Wrapper function for MAFFT multiple sequence alignment program
    
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
    return aquery, aref


def pair_align(query, ref, binpath='mafft'):
    """ Align query sequence to reference and drop insertions """
    aquery, aref = mafft(query, ref, binpath=binpath)
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
        "outfile", type=str,
        help="output, path to write aligned FASTA; if running MPI, "
             "script will append an integer suffix, e.g., '.1'"
    )
    parser.add_argument(
        "-b", "--binpath", default="mafft",
        help="Path to MAFFT executable file.  Defaults to 'mafft'."
    )
    parser.add_argument(
        "-f", "--format", default="fasta",
        help="format of sequence file format, passed to SeqIO.parse"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="print progress monitoring messages to stderr"
    )
    args = parser.parse_args()

    # read reference sequence
    reference = SeqIO.read(args.ref, args.format)
    ref = str(reference.seq)

    records = SeqIO.parse(args.infile, args.format)
    if nprocs == 1:
        outfile = open(args.outfile, 'w')
        for record in records:
            if args.verbose:
                sys.stderr.write(record.name+'\n')
                sys.stderr.flush()
            query = str(record.seq)
            aligned = mafft(query=query, ref=ref, binpath=args.binpath)
            outfile.write(f">{record.description}\n{aligned}\n")
    elif nprocs > 1:
        outfile = open(f"{args.outfile}.{my_rank}", 'w')
        for rn, record in enumerate(records):
            if rn % nprocs != my_rank:
                continue
            if args.verbose:
                sys.stderr.write(f"({my_rank}/{nprocs}) running {record.name}\n")
                sys.stderr.flush()
            aligned = mafft(query=str(record.seq), ref=ref, binpath=args.binpath)
            outfile.write(f">{record.description}\n{aligned}\n")

    outfile.close()
