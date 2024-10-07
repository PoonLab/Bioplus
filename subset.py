import random
from Bio import SeqIO
import argparse
import sys


description = """
Generate a random sample of a sequence dataset in a memory-
efficient way.  If neither -p nor -n are specified, prints
the number of sequences and exits.
"""


def sampler(infile, format='fasta', prop=0, num=0):
    """
    Iterate twice through file, first to count the number of sequence 
    records, and second to generate a random sample of records.  If both 
    num and prop are set to 0, print number of sequences to stderr and 
    return.

    :param infile:  str, path to file containing sequences
    :param format:  str, valid file format for Bio.SeqIO
    :param prop:  float, proportion of sequences to sample.  Used if
                  num = 0.
    :param num:  int, number of sequences to sample.

    :return:  None or list of SeqRecords
    """
    records = SeqIO.parse(infile, format)
    count = 0
    for _ in records:
        count += 1

    if num == 0:
        if prop == 0:
            sys.stdout.write(f"{count}\n")
            return None
        else:
            num = int(round(count * prop))
    
    idx = random.sample(range(count), num)
    records = SeqIO.parse(infile, format)
    sample = [record for i, record in enumerate(records) if i in idx]

    return sample


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description)
    parser.add_argument("infile", type=str,
                        help="Input, path to file containing sequences.")
    parser.add_argument("-f", "--format", type=str, default='fasta',
                        help="Option, format of input file (default: 'fasta')")
    parser.add_argument("-p", "--prop", type=float, default=0.,
                        help="Proportion of sequences to sample.")
    parser.add_argument("-n", "--num", type=int, default=0,
                        help="Number of sequences to sample.")
    parser.add_argument("-o", "--output", type=argparse.FileType('w'),
                        default=sys.stdout, 
                        help="Output, path to write down-sampled sequences "
                             "(default: stdout).")
    args = parser.parse_args()

    if args.prop < 0 or args.prop > 1:
        sys.stderr.write("ERROR: -p/--prop must be in interval [0,1]!\n")
        sys.exit()
    if args.num < 0:
        sys.stderr.write("ERROR: -n/--num cannot be negative!\n")
        sys.exit()

    sample = sampler(args.infile, args.format, prop=args.prop, num=args.num)
    if sample:
        SeqIO.write(sample, args.output, args.format)
