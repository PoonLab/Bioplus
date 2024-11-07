from Bio import SeqIO
import argparse
import sys


description = """
Generate a consensus sequence from a set of aligned sequences.
While this is possible in BioPython, it is rather convoluted.
"""

mixture_dict = {'W': 'AT', 'R': 'AG', 'K': 'GT', 'Y': 'CT', 'S': 'CG',
                'M': 'AC', 'V': 'AGC', 'H': 'ATC', 'D': 'ATG',
                'B': 'TGC', 'N': 'ATGC', '-': 'ATGC'}
ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.items())


def get_columns(aln, is_nuc=False):
    """ 
    Count frequencies of characters in each position 
    :param aln:  Bio.AlignIO object
    :param is_nuc:  bool, resolve nucleotide mixtures
    :return:  list of dicts (char, count) for all alignment positions
    """
    columns = None
    seqlen = None
    for record in aln:
        # handle first record
        if not columns:
            columns = [{nt: 1} for nt in record.seq]
            seqlen = len(record.seq)
            continue
        
        if len(record.seq) != seqlen:
            sys.stderr.write("ERROR: sequences are not the same length\n")
            sys.exit()
        
        for i, nt in enumerate(record.seq):
            if is_nuc and nt in mixture_dict:
                # count fractional stateseach
                chars = mixture_dict[nt]
                for ch in chars:
                    if ch not in columns[i]:
                        columns[i].update({ch: 0})
                    columns[i][ch] += 1./len(chars)
            else:            
                if nt not in columns[i]:
                    columns[i].update({nt: 0})
                columns[i][nt] += 1
    return columns


def conseq(columns, thresh=0.5, is_nuc=False):
    """
    Generate consensus sequence from an alignment based on 
    column-specific frequencies.
    
    :param columns:  list of dicts from get_columns()
    :param thresh:  float, minimum frequency for character to be included in mixture
    :param is_nuc:  bool, if True then handle multiple states as a nucleotide mixture
    :return:  str, plurality consensus sequence
    """
    seq = ''
    columns = get_columns(aln, is_nuc=is_nuc)
    for col in columns:
        total = sum(col.values())
        intermed = [(count/total, nt) for nt, count in col.items()]
        intermed.sort(reverse=True)
        states = []
        for freq, nt in intermed:
            if freq < thresh:
                break
            states.append(nt)
            
        states = ''.join(sorted(states))
        if len(states) == 1:
            seq += states
        else:
            if is_nuc and states in ambig_dict:
                seq += ambig_dict[states]  # convert nucleotide mixture
            else:
                if states[0] == '-':
                    seq += '-'
                else:
                    seq += '?'
    return seq
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="input, path to file containing aligned sequences")
    parser.add_argument("-f", "--format", default='fasta',
                        help="option, alignment file format (default: 'fasta')")
    parser.add_argument("-n", "--nuc", action="store_true",
                        help="Use nucleotide mixtures")
    parser.add_argument("-t", "--thresh", type=float, default=0.5,
                        help="If minority character is above this value (default 0.5), "
                        "then report '?' or mixture if '-n' is set.")
    args = parser.parse_args()
    
    aln = SeqIO.parse(args.infile, format=args.format)
    cseq = conseq(aln, thresh=args.thresh, is_nuc=args.nuc)
    sys.stdout.write(f"{cseq}\n")
