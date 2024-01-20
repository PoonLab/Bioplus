import argparse
from Bio import AlignIO, SeqIO, Phylo
import sys

description = """\
Down-sample a sequence alignment by building a tree and progressively
removing the shortest tips until only a target number of tips remain.
The tip labels are used to select sequences from the alignment. 
"""

def prunetree(phy, target):
    """
    Progressively remove the shortest terminal branches in the tree until
    we reach a target number of tips.
    :param phy:  Bio.Phylo object
    :param target:  int, number of tips we want to prune down to.
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
        "-n", "--target", type=int, required=True,
        help="Target number of sequences to reduce alignment to."
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

    # perform pruning and write resulting sequences to file
    pruned = prunetree(phy, target=args.target)
    for tip in pruned.get_terminals():
        record = records[tip.name]
        _ = SeqIO.write(record, args.outfile, "fasta")

