from Bio import Entrez, SeqIO
import sys
import time
import argparse
import re
from csv import DictWriter

description = """\
Retrieve metadata (e.g., sample collection dates) associated with sequences 
in the input file, based on their respective NCBI Genbank accession numbers.
"""

def find_accn(label, pat=re.compile("[A-Z]{1,3}[0-9]{5,8}\\.[0-9]|[NXAYW][CGTW"
                                    "ZMRP]_[A-Z0-9]+\\.[0-9]")):
    """
    Extract Genbank accession number from an input string.
    :param label:  str, string to process
    :param pat:  re.Pattern object
    :return:  first matching substring, or None if no matches
    """
    hits = pat.findall(label)
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
        help="Path to file containing sequences. Labels MUST contain Genbank "
             "accession numbers."
    )
    parser.add_argument(
        "--email", type=str, required=True,
        help="E-mail address for NCBI transactions."
    )
    parser.add_argument(
        "-f", "--format", type=str, default='fasta',
        help="Format of input sequence file; must be supported by Bio.SeqIO.  "
             "Defaults to 'fasta'.")
    parser.add_argument(
        "--db", type=str, default="nucleotide",
        help="Database to query, default 'nucleotide'."
    )
    parser.add_argument(
        "--delay", type=float, default=3.0,
        help="Number of seconds to pause between queries."
    )
    parser.add_argument(
        "--batch", type=int, default=50,
        help="Number of records to retrieve per request."
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType('w'), default=sys.stdout,
        help="option, path to write CSV output. Default is stdout."
    )
    args = parser.parse_args()

    assert '@' in args.email, "Use a valid e-mail address!"
    assert args.delay >= 1.0, "Delay must be at least 1.0 seconds!"
    assert args.batch > 0, "Batch size must be positive!"

    fields = ('organism', 'mol_type', 'isolate', 'isolation_source', 'host',
              'country', 'collection_date')
    Entrez.email = args.email

    # extract accessions from input
    records = SeqIO.parse(args.infile, args.format)
    headers = [record.description for record in records]
    intermed = [find_accn(header) for header in headers]

    accns = [a for a in intermed if a is not None]
    if len(accns) == 0:
        sys.stderr.write("\nERROR: Failed to parse any accession numbers from "
                         f"{args.infile.name}!\n\n")
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
    writer = DictWriter(args.outfile, fieldnames=('header',)+fields)
    writer.writeheader()

    # retrieve the metadata and write to file as CSV
    for i in range(0, len(accns), args.batch):
        query = ','.join(accns[i:(i+args.batch)])
        for row in get_metadata(query, fields, db=args.db):
            row.update({"header": headers[i]})
            writer.writerow(row)
        time.sleep(args.delay)
