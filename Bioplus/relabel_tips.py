description = """
Given a tree where some tips have missing metadata (e.g., subtype, country), 
impute these values by propagating labels from neighbouring tips.
"""

import argparse
import sys
import random
from Bio import Phylo


def climb(tips, curnode, pathlen, cutoff):
    """ Recursive function for traversing up a tree """
    pathlen += curnode.branch_length
    if pathlen < cutoff:
        if curnode.is_terminal():
            tips.append( (pathlen, curnode.name) )
        else:
            for c in curnode.clades:
                tips = climb(tips, c, pathlen, cutoff)
    return tips


def nearest(curnode, cutoff, parents):
    """
    Find nearest tips to a given tip in the tree.
    
    :param curnode:  Bio.Phylo.Clade object to search from
    :param cutoff:  float, maximum search distance in the tree
    :param parents:  dict, cached parent nodes
    
    :return: list, references to Clade objects
    """
    tips = []
    pathlen = curnode.branch_length
    p = parents[curnode]

    while p in parents:
        if pathlen >= cutoff:
            break
        for c in p.clades:
            if c == curnode:
                continue  # don't double back
            if c.is_terminal():
                this_dist = pathlen+c.branch_length
                if this_dist < cutoff:
                    tips.append( (this_dist, c.name) )
            else:
                tips.extend(climb([], c, pathlen, cutoff))
        
        curnode = p
        pathlen += p.branch_length
        p = parents[curnode]
    
    return tips


def relabel_tips(nwkfile, cutoff, delim, field, match='', 
                 ignore_case=False, k=1, verbose=False, **kwargs):
    """
    Use k-nearest neighbours in the tree to impute the missing values.
    If k=1 then use the closest non-missing tip label.  Otherwise use consensus.
    
    :param nwkfile:  _io.TextIOWrapper, handle (read-only) to file with Newick
                    tree string
    :param delim:  str, separates values in tip label
    :param field:  int, index of value to process
    :param match:  str, string to find exact match
    :param ignore_case:  bool, if True then convert both strings to lowercase
    :param k:  int, number of nearest neighbours to propagate labels from (default: 1)
    :param verbose:  bool, if True then print debugging messages
    
    :return:  Bio.Phylo.BaseTree object with imputed labels
    """
    ic = (str.lower if ignore_case else str)  # conditional function
    match = ic(match)
    
    phy = Phylo.read(nwkfile, 'newick')  # one tree only!
    
    # store parents of all nodes
    parents = {}
    for node in phy.find_clades(order='level'):
        for child in node.clades:
            parents.update({child: node})
    
    # locate all tips with missing labels and cache non-missing
    missing = set()
    labels = dict()
    tips = phy.get_terminals()
    for tip in tips:
        values = tip.name.split(delim)
        try:
            val = values[field]
        except IndexError:
            sys.stderr.write(f"Failed to locate field {field} in header: {values}\n")
            continue
        
        if ic(val) == match:
            missing.add(tip)
        else:
            labels.update({tip.name: val})
    
    if verbose:
        sys.stderr.write(f"Found {len(missing)} tips with missing labels.\n")
    
    for tip in missing:
        neighbours = nearest(tip, cutoff, parents)
        neighbours = [(l, n) for l, n in neighbours if n in labels]
        neighbours.sort()  # order of increasing path length
        
        knn = neighbours[:k]
        if len(knn) == 0:
            sys.stderr.write(f"\nERROR: no neighbours for {tip.name}!\n"
                             "Try increasing --cutoff.\n\n")
            sys.exit()
        
        if k > 1:
            poll = [labels[n] for _, n in knn]
            newval = consensus(poll)
        else:
            try:
                newval = labels[knn[0][1]]  # I'm feeling lucky
            except IndexError:
                sys.stderr.write(f"{knn}\n")
                raise

        # update tip label
        values = tip.name.split(delim)
        old_name = tip.name
        tip.name = delim.join(values[:field]+[newval]+values[(field+1):])
        
        if verbose:
            sys.stderr.write(f"{old_name} -> {tip.name}\n")

    return phy


def consensus(poll):
    """ Determine the consensus from a poll of labels from neighbours """
    intermed = [(poll.count(l), l) for l in set(poll)]
    intermed.sort(reverse=True)
    counts = [c for c, l in intermed]
    try:
        max_count = counts[0]
    except IndexError:
        sys.stderr.write(f"poll: {poll}\n")
        raise
    
    if counts.count(max_count) > 1:
        vals = [l for c, l in intermed if c == max_count]
        return random.sample(vals, 1)  # random resolution of tie
    else:
        return intermed[0][1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "nwkfile", 
        type=argparse.FileType('r'), 
        help="Path to file containing Newick tree.")
    parser.add_argument(
        '-o', "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="Path to write output file (default: stdout).")
    parser.add_argument(
        "-d", "--delim", type=str, default="_",
        help="Delimiting character that separates values in tip label \
            (default: '_').")
    parser.add_argument(
        "-f", "--field", type=int, default=-1,
        help="Integer (0-)index to field containing target metadata \
            (default: -1).")
    parser.add_argument(
        "-m", "--match", type=str, default="",
        help="String that identifies missing values (default: '').")
    parser.add_argument(
        '-i', "--ignore_case", action="store_true", default=False,
        help="Ignore case when matching search string.")
    parser.add_argument(
        '--cutoff', type=float, default=0.1,
        help="Search distance in units of branch length (default: 0.1).")
    parser.add_argument(
        '-k', type=int, default=1, 
        help="Number of nearest neighbours to impute labels from.")
    parser.add_argument(
        "--format", default="newick", 
        help="File format for Phylo.read (default: 'newick')")
    parser.add_argument(
        "--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()
    
    tree = relabel_tips(**vars(args))  # unpack Namespace to keyword arguments
    Phylo.write(tree, args.outfile, format='newick')
