import argparse
from Bio import AlignIO, SeqIO, Phylo
import sys

description = """\
Down-sample a sequence alignment by building a tree and progressively
removing the shortest tips until only a target number of tips remain.
The tip labels are used to select sequences from the alignment. 
"""


def prune_tips(phy, target, cache=False):
    """
    Progressively remove the shortest terminal branches in the tree until
    we reach a target number of tips.
    :param phy:  Bio.Phylo object
    :param target:  float, number of tips we want to prune down to.
    :param cache:  if True, copy label of pruned tip to closest terminal node
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


def prune_length(phy, target):
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
        _ = phy.prune(tips[0])
        tlen -= tips[0].branch_length
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
        "-n", "--target", type=float, required=True,
        help="Target number of sequences (default) or tree length "
             "(--length) to reduce alignment to."
    )

    parser.add_argument(
        "--seq", type=argparse.FileType('r'), default=None,
        help="Path to file containing sequence alignment.  Script will "
             "output reduced set of sequences instead of tree."
    )

    parser.add_argument(
        "--length", action="store_true",
        help="Prune tree to target total length instead of tip count."
    )
    parser.add_argument(
        "-f", "--format", default="fasta", type=str,
        help="Specify format of input alignment.  Must be supported by "
             "Bio.AlignIO (default 'fasta').  Only used for --seq."
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType('w'),
        help="Path to write down-sampled alignment FASTA (--seq) or tree."
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
        aln = AlignIO.read(args.infile, args.format)
        records = dict([(record.name, record) for record in aln])
        labels = set(records.keys())
        if tip_names != labels:
            sys.stderr.write("ERROR: Input tree labels do not match alignment.")
            sys.exit()

    # perform pruning
    if args.length:
        tlen = phy.total_branch_length()
        sys.stderr.write(f"Starting tree length: {tlen}\n")
        if args.target <= 0:
            sys.stderr.write("ERROR! Target length must be positive!\n")
            sys.exit()
        pruned = prune_length(phy, target=args.target)
    else:
        sys.stderr.write(f"Starting tip count: {len(phy.get_terminals())}\n")
        pruned = prune_tips(phy, target=args.target,
                            cache=args.csvfile is not None)

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
            for tip in pruned.get_terminals():
                if not hasattr(tip, 'cache'):
                    continue
                for label in tip.cache:
                    args.csvfile.write(f"{tip.name},{label}\n")
