import argparse
from Bio import AlignIO, SeqIO, Phylo
import sys
import csv

description = """\
Down-sample a sequence alignment by building a tree and progressively
removing the shortest tips until only a target number of tips remain.
The tip labels are used to select sequences from the alignment. 
"""


def prune_tips(phy, target, cache=False, trace=False):
    """
    Progressively remove the shortest terminal branches in the tree until
    we reach a target number of tips.
    :param phy:  Bio.Phylo object
    :param target:  float, number of tips we want to prune down to.
    :param cache:  bool, copy label of pruned tip to closest terminal node
    :param trace:  bool, print tree length as a function of number of tips
    :return:  Bio.Phylo object
    """
    tips = phy.get_terminals()  # returns a list of Clade objects
    if target >= len(tips):
        sys.stderr.write(f"prunetree: requirement already met "
                         f"({len(tips)}<={target})\n")
        return phy

    while len(tips) > target:
        # find shortest tip
        # FIXME: faster to update instead of resorting every time?
        tips = sorted(tips, key=lambda x: x.branch_length)
        if trace:
            sys.stdout.write(f"{len(tips)},{phy.total_branch_length()}\n")
            sys.stdout.flush()
        tip = tips[0]
        parent = phy.prune(tip)
        if cache:
            kin = parent.get_terminals(order="level")[0]
            if not hasattr(kin, "cache"):
                kin.cache = []
            kin.cache.append(tip.name)
            if hasattr(tip, "cache"):
                kin.cache.extend(tip.cache)

        tips = tips[1:]  # instead of calling get_terminals() again

    return phy


def prune_length(phy, target, cache=False):
    """
    Progressively remove the shortest terminal branches in the tree until
    we reach a target tree length.
    :param phy:  Bio.Phylo object
    :param target:  float, target tree length to prune to
    :return:  Bio.Phylo object
    """
    tlen = phy.total_branch_length()
    if target >= tlen:
        sys.stderr.write(f"prune_length: requirement already met "
                         f"({tlen}<={target}")
    tips = phy.get_terminals()
    while tlen > target:
        # we have to re-sort every time because removing a branch
        # will lengthen another branch
        tips = sorted(tips, key=lambda x: x.branch_length)
        tip = tips[0]
        parent = phy.prune(tip)
        if cache:
            kin = parent.get_terminals(order="level")[0]
            if not hasattr(kin, "cache"):
                kin.cache = []
            kin.cache.append(tip.name)
            if hasattr(tip, "cache"):
                kin.cache.extend(tip.cache)
        tlen -= tip.branch_length
        tips = tips[1:]  # update list
    return phy


def prune_tiplen(phy, target, cache=False):
    """
    Progressively remove terminal branches that have a length below
    the target minimum length (--mode tiplen).
    :param phy:  Bio.Phylo object
    :param target:  float, minimum tip length
    :return:  Bio.Phylo object
    """
    tips = phy.get_terminals()
    while True:
        # we have to re-sort every time because removing a branch
        # will lengthen another branch
        tips = sorted(tips, key=lambda x: x.branch_length)
        tip = tips[0]
        if tip.branch_length > target:
            break

        parent = phy.prune(tip)
        if cache:
            kin = parent.get_terminals(order="level")[0]
            if not hasattr(kin, "cache"):
                kin.cache = []
            kin.cache.append(tip.name)
            if hasattr(tip, "cache"):
                kin.cache.extend(tip.cache)

        tips = tips[1:]  # update list
    return phy


if __name__ == "__main__":
    # command-line interface
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "tree", type=argparse.FileType('r'),
        help="Path to file containing Newick tree string."
    )
    parser.add_argument(
        "-t", "--target", type=float, default=None,
        help="Target quantity, conditional on --mode; e.g., -t 100 "
             "--mode ntips will prune the tree down to 100 tips. "
             "If this is not specified, then script will provide info "
             "about input tree."
    )
    parser.add_argument(
        "--seq", type=argparse.FileType('r'), default=None,
        help="Path to file containing sequence alignment.  Script will "
             "output reduced set of sequences instead of tree."
    )
    parser.add_argument(
        "--mode", choices=['ntips', 'treelen', 'tiplen'], default='ntips',
        help="Prune tree to a target of [ntips] number of tips; "
             "[treelen] total branch length of tree; or [tiplen] "
             "minimum branch length."
    )
    parser.add_argument(
        "-f", "--format", default="fasta", type=str,
        help="Specify format of input alignment.  Must be supported by "
             "Bio.AlignIO (default 'fasta').  Only used for --seq."
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="Path to write down-sampled alignment FASTA (--seq) or tree."
             "Defaults to stdout."
    )
    parser.add_argument(
        "--csvfile", type=argparse.FileType('w'), default=None,
        help="Optional, record the labels of pruned tips associated with "
             "remaining tips into a CSV file."
    )
    args = parser.parse_args()

    phy = Phylo.read(args.tree, "newick")
    tip_names = set([tip.name for tip in phy.get_terminals()])

    if args.seq:
        # checks whether sequences have the same length
        aln = AlignIO.read(args.seq, args.format)
        records = dict([(record.description.split("'")[0], record) for record in aln])
        labels = set(records.keys())
        if tip_names != labels:
            sys.stderr.write("ERROR: Input tree labels do not match alignment.\n")
            sys.stderr.write(f"{tip_names.difference(labels)}\n")
            sys.exit()

    # perform pruning
    if args.mode == "treelen":
        tlen = phy.total_branch_length()
        if args.target is None or args.target <= 0:
            sys.stderr.write(f"Starting tree length: {tlen}\n")
            sys.exit()
        pruned = prune_length(phy, target=args.target, cache=args.csvfile is not None)
    elif args.mode == "tiplen":
        if args.target is None:
            tips = phy.get_terminals()
            tiplens = sorted([tip.branch_length for tip in tips])
            sys.stderr.write(f"Shortest tip lengths:\n")
            sys.stderr.write(f"  min:    {tiplens[0]}\n")
            sys.stderr.write(f"  2.5%:   {tiplens[round(0.025 * len(tips))]}\n")
            sys.stderr.write(f"  25%:    {tiplens[round(0.25 * len(tips))]}\n")
            sys.stderr.write(f"  median: {tiplens[round(0.5*len(tips))]}\n")
            sys.exit()
        pruned = prune_tiplen(phy, target=args.target, cache=args.csvfile is not None)
    elif args.mode == "ntips":
        if args.target is None:
            prune_tips(phy, target=3, trace=True)
            #sys.stderr.write(f"Starting tip count: {len(phy.get_terminals())}\n")
            sys.exit()
        pruned = prune_tips(phy, target=args.target,
                            cache=args.csvfile is not None)
    else:
        sys.stderr.write(f"ERROR: Unknown --mode option {args.mode}, exiting\n")
        sys.exit()

    if args.seq:
        # write resulting sequences to output file
        for tip in pruned.get_terminals():
            record = records[tip.name]
            _ = SeqIO.write(record, args.outfile, "fasta")
    else:
        # write pruned tree to output file
        Phylo.write(pruned, args.outfile, 'newick')

    if args.csvfile:
        # write pruned labels to CSV
        writer = csv.writer(args.csvfile)
        for tip in pruned.get_terminals():
            if not hasattr(tip, 'cache'):
                continue
            for label in tip.cache:
                writer.writerow([tip.name, label])
