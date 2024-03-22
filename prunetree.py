import argparse
from Bio import AlignIO, SeqIO, Phylo
import sys

description = """\
Down-sample a sequence alignment by building a tree and progressively
removing the shortest tips until only a target number of tips remain.
The tip labels are used to select sequences from the alignment. 
"""

def prune_tips(phy, target):
    """
    Progressively remove the shortest terminal branches in the tree until
    we reach a target number of tips.
    :param phy:  Bio.Phylo object
    :param target:  float, number of tips we want to prune down to.
    :return:  Bio.Phylo object
    """
    tips = phy.get_terminals()  # returns a list of Clade objects
    if target >= len(tips):
        sys.stderr.write(f"prunetree: requirement already met "
                         f"({len(tips)}<={target})\n")
        return phy

    while len(tips) > target:
        # find shortest tip
        # FIXME: faster to update instead of rebuilding every time?
        tips = sorted(tips, key=lambda x: x.branch_length)
        _ = phy.prune(tips[0])
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
        "infile", type=argparse.FileType('r'),
        help="Path to file containing sequence alignment."
    )
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
        "--length", action="store_true",
        help="Prune tree to target total length instead of tip count."
    )
    parser.add_argument(
        "-f", "--format", default="fasta", type=str,
        help="Specify format of input alignment.  Must be supported by "
             "Bio.AlignIO (default 'fasta')."
    )
    parser.add_argument(
        "-o", "--outfile", default=sys.stdout, type=argparse.FileType('w'),
        help="Path to write down-sampled alignment FASTA file."
    )
    args = parser.parse_args()

    # checks whether sequences have the same length
    aln = AlignIO.read(args.infile, args.format)
    records = dict([(record.name, record) for record in aln])
    labels = set(records.keys())

    # check whether tree matches alignment
    phy = Phylo.read(args.tree, "newick")
    tip_names = set([tip.name for tip in phy.get_terminals()])
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
        pruned = prune_tips(phy, target=args.target)

    # write resulting sequences to output file
    for tip in pruned.get_terminals():
        record = records[tip.name]
        _ = SeqIO.write(record, args.outfile, "fasta")

