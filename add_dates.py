import argparse
from Bio import SeqIO
from csv import DictReader
from datetime import date
import re


description = """\
Use metadata to relabel sequences in the associated FASTA file.
Optionally, use the dates to down-sample the sequences to a target 
size so that dates are evenly distributed.
"""

def lubridate(string, delim=',|-|/| '):
    """
    Guess the date format of a string, starting with ISO format.
    :param string:  str, serialization of a calendar date in unknown format
    :param delim:  str, a regular expression comprising all potential characters
                   delimiting date fields, separated by "|"
    :return:  datetime.date object
    """
    try:
        # first try ISO format
        return date.fromisoformat(string)
    except ValueError:
        items = re.split(delim, string)
        if len(items) == 1:
            try:
                return int(string)  # year?
            except ValueError:
                return None  # failed
        elif len(items) == 2:
            # try some month-year combos
            
        elif len(items) == 3:

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="File of sequences to process")
    parser.add_argument("metadata", type=argparse.FileType('r'),
                        help="CSV file containing metadata (from get_metadata.py)")
    parser.add_argument("--field", type=str, default="collection_date",
                        help="Column label for field in CSV containing dates.")
    parser.add_argument("-f", "--format", type=str, default="fasta",
                        help="Format specifier for input sequence file.")
    args = parser.parse_args()

    # load metadata and parse collection dates
    coldates = {}
    for row in DictReader(args.metadata):
        coldates.update({row["header"]: row[args.field]})

    headers, seqs = [], []
    for record in SeqIO.parse(args.infile, args.format):
        headers.append(record.description)
        seqs.append(record)
