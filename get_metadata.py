from Bio import Entrez, SeqIO, Phylo
import sys
import time
import argparse
import re
from csv import DictWriter

description = """\
Retrieve metadata (e.g., sample collection dates) associated with sequences 
in the input file, based on their respective NCBI Genbank accession numbers.
"""

pat1 = re.compile("[A-Z]{1,3}[0-9]{5,8}\\.[0-9]|[NXAYW][CGTWZMRP]_[A-Z0-9]+\\.[0-9]")
pat2 = re.compile("[A-Z]{1,3}[0-9]{5,8}\\.[0-9]|[NXAYW][CGTWZMRP]_[A-Z0-9]+")


def find_accn(label, versioned=True):
    """
    Extract Genbank accession number from an input string.
    :param label:  str, string to process
    :param versioned:  bool, if True then accession number must have ".[0-9]" suffix
    :return:  first matching substring, or None if no matches
    """
    if versioned:
        hits = pat1.findall(label)
    else:
        hits = pat2.findall(label)
    if len(hits) == 0:
        return None
    return hits[0]


def get_metadata(accn, fields, db='nucleotide'):
    """
    Retrieve source metadata associated with sequence for a given accession
    number.
    :param accn:  str, Genbank accession number.  Multiple numbers may be
                  submitted as a comma-delimited string.
    :param fields:  tuple, field names to extract from source qualifier, e.g.,
                    ("country", "collection_date").
    :param db:  str, Genbank database to query
    :return:  dict, values associated with fields
    """
    handle = Entrez.efetch(db=db, rettype='gb', retmode='text', id=accn)
    output = SeqIO.parse(handle, format='genbank')
    for record in output:
        source = record.features[0]
        values = [source.qualifiers.get(f, [''])[0] for f in fields]
        yield dict(zip(fields, values))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "infile", type=argparse.FileType('r'),
        help="Path to file containing sequences or tree. Labels MUST contain Genbank "
             "accession numbers."
    )
    parser.add_argument(
        "--email", type=str, required=True,
        help="E-mail address for NCBI transactions."
    )
    parser.add_argument(
        "-f", "--format", type=str, default='fasta',
        help="Format of input sequence/tree file; must be supported by Bio.SeqIO or "
             "Phylo.  Defaults to 'fasta' (use 'newick' for trees).")
    parser.add_argument(
        "--db", type=str, default="nucleotide",
        help="Database to query, default 'nucleotide'."
    )
    parser.add_argument(
        "--delay", type=float, default=3.0,
        help="Number of seconds to pause between queries (default 3)."
    )
    parser.add_argument(
        "--batch", type=int, default=50,
        help="Number of records to retrieve per request (default 50)."
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="option, path to write CSV output.  If --relabel, then "
             "outputs a new sequence or tree file.  Default is stdout."
    )
    parser.add_argument(
        "--nover", action="store_true",
        help="Set if accession number has no version number suffix, e.g., "
             "AB123456 instead of AB123456.1"
    )
    parser.add_argument(
        "--relabel", action="store_true",
        help="Use metadata to replace sequence/tree labels with: "
             "[organism]_[host]_[country]_[collection date]."
    )
    args = parser.parse_args()

    assert '@' in args.email, "Use a valid e-mail address!"
    assert args.delay >= 1.0, "Delay must be at least 1.0 seconds!"
    assert args.batch > 0, "Batch size must be positive!"

    fields = ('organism', 'mol_type', 'isolate', 'isolation_source', 'host',
              'geo_loc_name', 'collection_date')
    Entrez.email = args.email

    # extract accessions from input
    if args.format == "newick":
        phy = Phylo.read(args.infile, "newick")
        headers = [tip.name for tip in phy.get_terminals()]
    else:
        records = SeqIO.parse(args.infile, args.format)
        headers = [record.description for record in records]

    intermed = [find_accn(header, versioned=not args.nover) for header in headers]
    accns = [a for a in intermed if a is not None]
    if len(accns) == 0:
        sys.stderr.write("\nERROR: Failed to parse any accession numbers from "
                         f"{args.infile.name}!\n")
        if args.infile.name.endswith(".nwk"):
            sys.stderr.write("*** This seems to be a tree, did you forget to set `--format "
                             "newick`? ***\n")
        sys.exit()

    sys.stderr.write(
        f"Extracted {len(accns)} accession numbers from input file:\n")
    sys.stderr.write(
        f"  {', '.join(accns[:3])}, ..., {', '.join(accns[-3:])}\n")
    if len(accns) < len(intermed):
        sys.stderr.write(
            f"Failed to parse accessions from {len(intermed)-len(accns)} "
            f"records.\n")
    sys.stderr.flush()

    # prepare output file
    if not args.relabel:
        writer = DictWriter(args.outfile, fieldnames=('header',)+fields)
        writer.writeheader()

    # retrieve the metadata and write to file as CSV
    labels = {}
    for i in range(0, len(accns), args.batch):
        query = ','.join(accns[i:(i+args.batch)])
        rows = get_metadata(query, fields, db=args.db)
        for j, row in enumerate(rows):
            if args.relabel:
                lab = f"{accns[i+j]}_{row['geo_loc_name']}_{row['collection_date']}"
                labels.update({headers[i+j]: lab})
            else:
                row.update({"header": headers[i+j]})
                writer.writerow(row)
        time.sleep(args.delay)

    if args.relabel:
        if args.format == 'newick':
            for tip in phy.get_terminals():
                tip.name = labels[tip.name]
            Phylo.write(phy, args.outfile, 'newick')
        else:
            # reload input file
            records = SeqIO.parse(args.infile, args.format)
            for record in records:
                record.description = labels[record.description]
                SeqIO.write(record, args.outfile, args.format)
    
    args.outfile.close()
