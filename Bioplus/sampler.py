from Bio import SeqIO
import argparse
import random
import sys
from datetime import datetime


def sample_by_date(records, target, delim="_", pos=-1, mode='year'):
    """ 
    Bin records by date intervals and draw random samples per bin
    up to a maximum number (max).
    @param records:  Bio.SeqIO iterable
    @param target:  int, maximum number of records to sample per time interval
    @param delim:  str, for parsing sequence labels
    @param pos:  int, position of date object in delimited sequence label
    @param mode:  str, binning by 'year' or by 'month'
    @return:  list, Bio.SeqIO SeqRecord objects
    """
    by_date = {}
    for record in records:
        date_str = record.description.split(delim)[pos]
        try:
            dt = datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError:
            try:
                dt = datetime.strptime(date_str, "%Y")
            except ValueError:
                sys.stderr.write(f"ERROR: Failed to parse date {date_str} - "
                                 f"skipping\n")
                continue
        
        if mode == 'year':
            key = dt.year
        elif mode == 'month':
            key = (dt.year, dt.month)
        else:
            sys.stderr.write(f"ERROR: Unrecognized mode '{mode}' in "
                             f"sample_by_date, exiting\n")
            sys.exit()
        
        if key not in by_date:
            by_date.update({key: []})
        by_date[key].append(record)

    sample = []
    for key, bin in by_date.items():
        if len(bin) <= target:
            sample.extend(bin)  # not enough to sample
        else:
            sample.extend(random.sample(bin, target))

    return sample


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Downsample sequences")
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="Path to file containing sequences to downsample.")
    parser.add_argument("target", type=int, help="Target number of sequences.")
    parser.add_argument("-f", "--format", type=str, default='fasta',
                        help="Input format, default 'fasta'")
    parser.add_argument("-o", "--output", type=str, default='fasta',
                        help="Output format, default 'fasta'")
    parser.add_argument("-m", "--method", type=str, default="random",
                        choices=["random", "head", "tail", "year", "month"],
                        help="Sampling method to use: 'random', sample "
                             "completely at random; 'head', write first n; "
                             "'tail', write last n; 'year', downsample by "
                             "year; 'month', downsample by month.")
    parser.add_argument("--delim", type=str, default="_",
                        help="Delimiter to split sequence labels for --method "
                             "'year' or 'month'")
    parser.add_argument("--pos", type=int, default=-1,
                        help="Position of date in sequence label, used only "
                             "for --method 'year' or 'month'.  Defaults to -1.")
    args = parser.parse_args()

    # load all sequences into memory
    records = list(SeqIO.parse(args.infile, args.format))
    if len(records) == 0:
        sys.stderr.write("Failed to read records, are you sure you set the correct "
              "--format?\n")
        sys.exit()
    sys.stderr.write(f"Read {len(records)} sequences from file\n")

    if args.method == "random":
        sample = random.sample(records, args.target)
    elif args.method == "head":
        sample = records[:args.target]
    elif args.method == "tail":
        sample = records[-args.target:]
    elif args.method in ["year", "month"]:
        sample = sample_by_date(records, target=args.target, delim=args.delim,
                                pos=args.pos, mode=args.method)
    else:
        print(f"Error: unrecognized method {args.method}")
        sys.exit()

    # write to stdout
    SeqIO.write(sample, sys.stdout, args.output)
