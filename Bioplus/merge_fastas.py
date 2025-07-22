from Bio import SeqIO
import argparse
import sys

description = """
Concatenate sequences in FASTA files specified by the user, in the order given
on the command line, by matching labels.  Raise a warning if matching labels 
are not found.  Write result to a new FASTA file (default: stdout).
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument("inputs", nargs='+', help="Any number of FASTA files to combine.")
parser.add_argument("-o", "--output", type=argparse.FileType('w'), 
                    default=sys.stdout, help="Output path (default: stdout)")
parser.add_argument("-f", "--force", action="store_true", 
                    help="Optional, output sequence even if no match found.\n")
args = parser.parse_args()

n = len(args.inputs)
if n < 2:
    sys.stderr.write("Less than two files given, exiting..\n")
    sys.exit()

seqdicts = []
for fp in args.inputs:
    records = SeqIO.parse(fp, 'fasta')
    seqdict = {}
    for record in records:
        seqdict.update({record.description: record.seq})
    seqdicts.append(seqdict)

for label, seq in seqdicts[0].items():
    failed = False
    for i in range(1, n):
        if label not in seqdicts[i]:
            failed = True
        else:
            seq += seqdicts[i][label]
    if failed and not args.force:
        continue
    args.output.write(f">{label}\n{seq}\n")
