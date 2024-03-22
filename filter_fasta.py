import argparse
from Bio import SeqIO
from csv import DictReader
import sys


description = """\
Filter sequences in a FASTA file based on the presence or 
absence of metadata.  Optionally append the metadata to 
sequence labels.
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="File of sequences to process")
    parser.add_argument("metadata", type=argparse.FileType('r'),
                        help="CSV file containing metadata (from get_metadata.py)")
    parser.add_argument("-f", "--format", type=str, default="fasta",
                        help="Format specifier for input sequence file.")
    parser.add_argument("--name", type=str, default="header",
                        help="Column label of labels to match metadata to FASTA. "
                             "Defaults to 'header'.")
    parser.add_argument("--field", type=str,
                        help="Column label of field in CSV to filter FASTA.")
    parser.add_argument("--append", action="store_true",
                        help="If set, add field to FASTA labels.")
    parser.add_argument("--sep", type=str, default="_",
                        help="Delimiter to update sequence label with metadata. "
                             "Used only with --append.")
    parser.add_argument("--missing", type=str, default="NA",
                        help="String that indicates a value is missing. "
                             "Defaults to 'NA'.")
    args = parser.parse_args()

    # parse CSV file
    values = {}
    reader = DictReader(args.metadata)
    for row in reader:
        values.update({row[args.name]: row[args.field]})

    # apply parsed dates to sequence labels
    for record in SeqIO.parse(args.infile, args.format):
        header = record.name
        value = values.get(header, None)
        if value is None:
            sys.stderr.write(f"Warning: Could not retrieve metadata for {header}, skipping.")
            continue

        if value == "" or value == args.missing:
            continue  # missing value, omit record

        # write record to output
        if args.append:
            res = f">{record.description}{args.sep}{value}\n{record.seq}\n"
        else:
            res = f">{record.description}\n{record.seq}\n"
        sys.stdout.write(res)
