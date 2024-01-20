from Bio import SeqIO
import argparse
import random
import sys

parser = argparse.ArgumentParser("Downsample sequences")
parser.add_argument("infile", type=argparse.FileType('r'),
                    help="Path to file containing sequences to downsample.")
parser.add_argument("target", type=int, help="Target number of sequences.")
parser.add_argument("-f", "--format", type=str, default='fasta',
                    help="Input format, default 'fasta'")
parser.add_argument("-o", "--output", type=str, default='fasta',
                    help="Output format, default 'fasta'")
parser.add_argument("-m", "--method", type=str, default="random",
                    choices=["random", "head", "tail"],
                    help="Sampling method to use: 'random', sample completely "
                         "at random; 'head', write first n; 'tail', write "
                         "last n; ")
args = parser.parse_args()

# load all sequences into memory
records = list(SeqIO.parse(args.infile, args.format))
if len(records) == 0:
    print("Failed to read records, are you sure you set the correct --format?")
    sys.exit()
print(f"Read {len(records)} sequences from file")

if args.method == "random":
    sample = random.sample(records, args.target)
elif args.method == "head":
    sample = records[:args.target]
elif args.method == "tail":
    sample = records[-args.target:]
else:
    print(f"Error: unrecognized method {args.method}")
    sys.exit()

# write to stdout
SeqIO.write(sample, sys.stdout, args.output)
