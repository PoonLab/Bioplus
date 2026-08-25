description = """
Given a tree where some tips have missing metadata (e.g., subtype, country), 
impute these values by propagating labels from neighbouring tips.
"""

import argparse
import sys
from Bio import Phylo


def climb(tips, curnode, pathlen, cutoff):
    """ Recursive function for traversing up a tree """
    pathlen += curnode.branch_length
    if pathlen < cutoff:
        if curnode.is_terminal():
            tips.append( (curnode.name, pathlen) )
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


def relabel_tips(nwkfile, cutoff, delim, field, match='', ignore_case=False, k=1):
    """
    Use k-nearest neighbours in the tree to impute the missing values.
    If k=1 then use the closest non-missing tip label.  Otherwise use consensus.
    
    :param tree:  object of class Bio.Phylo.BaseTree
    :param delim:  str, separates values in tip label
    :param field:  int, index of value to process
    :param match:  str, string to find exact match
    :param ignore_case:  bool, if True then convert both strings to lowercase
    
    :return:  Bio.Phylo.BaseTree object with imputed labels
    """
    ic = (str.lower if ignore_case else str)  # conditional function
    match = ic(match)
        
    phy = Phylo.read(nwkfile, 'newick')
    
    # store parents of all nodes
    parents = {}
    for node in phy.find_clades(order='level'):
        for child in node.clades:
            parents.update({child: node})
    
    # locate all tips with missing labels
    missing = set()
    tips = phy.get_terminals()
    for tip in tips:
        values = tip.name.split(delim)
        try:
            val = ic(values[field])
        except IndexError:
            sys.stderr.write(f"Failed to locate field {field} in header: {values}\n")
            continue
        
        if val == match:
            missing.add(tip)
    
        
    # otherwise locate nearest neighbours
    neighbours = nearest(tip, cutoff, parents)
    neighbours.sort()  # order of increasing path length
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "nwkfile", type=argparse.FileType('r'), help="Path to file containing Newick tree.")
    parser.add_argument(
        '-o', "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="Path to write output file (default: stdout)."
        )
    parser.add_argument(
        "-d", "--delim", type=str, default="_",
        help="Delimiting character that separates values in tip label (default: '_')."
        )
    parser.add_argument(
        "-f", "--field", type=int, default=-1,
        help="Integer (0-)index to field containing target metadata (default: -1)."
    )
    parser.add_argument(
        "-p", "--pattern", type=str, default="",
        help="String that identifies missing values (default: '')."
    )
    parser.add_argument(
        '-i', "--ignorecase", action="store_true", default=False,
        help="Ignore case when matching search string.")
    parser.add_argument(
        '--cutoff', type=float, default=0.1,
        help="Search distance in units of branch length (default: 0.1)."
    )
    args = parser.parse_args()
    
    relabel_tips(nwkfile=args.nwkfile, delim=args.delim, field=args.field, 
                 match=args.pattern, ignore_case=args.ignorecase, cutoff=args.cutoff)
    