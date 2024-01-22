# Bioplus
### A collection of Python scripts extending the functionality of Biopython

* `add_dates.py` - IN PROGRESS - Use metadata to relabel sequences in the associated FASTA file. Optionally, use the dates to down-sample the sequences to a target size so that dates are evenly distributed.
* `dist.py` - Calculate a k-mer distance (Euclidean, spectrum kernel), p-distance or Jukes-Cantor corrected distance for an alignment of nucleotide sequences.
* `get_metadata.py` - Retrieve metadata (e.g., sample collection dates) associated with sequences in the input file, based on their respective NCBI Genbank accession numbers.
* `permute_fasta.py` - Random permutation of FASTA file. This script reads a FASTA-formatted file containing aligned sequences and applies a random permutation of alignment columns (nucleotide sites), writing the result to a new FASTA file.  This has the benefit of providing some data security (by scrambling the original sequences) while preserving information about evolutionary relationships.
* `prunetree.py` - Down-sample a sequence alignment by building a tree and progressively removing the shortest tips until only a target number of tips remain. The tip labels are used to select sequences from the alignment.
* `sampler.py` - Reduce the number of sequences in a file by random sampling.
