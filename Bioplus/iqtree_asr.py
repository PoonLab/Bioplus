"""
Process a sequence alignment and ROOTED tree for ancestral sequence
reconstruction using IQ-TREE, preserving the root node.
"""

import argparse
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import subprocess
import tempfile
from os.path import basename, splitext


parser = argparse.ArgumentParser()
parser.add_argument("aln", type=argparse.FileType('r'),
                    help="<input> file containing alignment")
parser.add_argument("nwk", type=argparse.FileType('r'),
                    help="<input> file containing Newick tree (rooted)")
parser.add_argument("--bin", default="iqtree2", type=str,
                    help="Path to IQ-TREE executable")
parser.add_argument("-m", "--model", default="HKY", type=str,
                    help="IQ-TREE substitution model (default HKY).")
parser.add_argument("--prefix", type=str, help="Output file prefix.")
parser.add_argument("-nt", type=int, help="Number of threads", default=1)
parser.add_argument("--redo", action="store_true", help="Overwrite previous outputs")
args = parser.parse_args()

if args.prefix is None:
    args.prefix = splitext(basename(args.nwk.name))[0]
    print(f"No prefix provided, defaulting to {args.prefix}")

# import tree
phy = Phylo.read(args.nwk, "newick")
if len(phy.root.clades) > 2:
    print("ERROR: This does not seem to be a rooted tree!")
    sys.exit()

# add a zero-length branch to root node
new_tip = Phylo.BaseTree.Clade(branch_length=0, name="ROOT")
phy.root.clades.append(new_tip)

# write to temp file
temptree = tempfile.NamedTemporaryFile(mode='w+', delete=False, encoding="utf-8")
Phylo.write(phy, temptree, "newick")
temptree.close()

# add empty sequence to alignment
aln = AlignIO.read(args.aln, "fasta")
empty_seq = Seq("-"*aln.get_alignment_length())
dummy = SeqRecord(empty_seq, id="ROOT", description="")
aln.append(dummy)

tempalign = tempfile.NamedTemporaryFile(mode="w+", delete=False, encoding="utf-8")
AlignIO.write(aln, tempalign, "fasta")
tempalign.close()

# call IQ-TREE
cmd = [args.bin, "-s", tempalign.name, "-te", temptree.name,
       "-m", args.model, "-asr", "--prefix", args.prefix, 
       "-nt", str(args.nt), "-keep_empty_seq"]
if args.redo:
    cmd.append("--redo")
subprocess.check_call(cmd)


