import argparse
from Bio import SeqIO
from csv import DictReader
from datetime import date, datetime
import re
import sys


description = """\
Use metadata to relabel sequences in the associated FASTA file.
Optionally, use the dates to down-sample the sequences to a target 
size so that dates are evenly distributed.  Output written to 
stdout stream.  "header" field of metadata CSV must match the 
sequence names in the FASTA file.
"""


def str2year(s):
    """ Convert string representation of year to integer value """
    if len(s) == 4:
        return datetime.strptime(s, "%Y").year
    elif len(s) == 2:
        return datetime.strptime(s, "%y").year
    else:
        return None


def str2month(s):
    """ Convert string representation of month to integer value"""
    if len(s) == 3:
        # abbreviated month
        return datetime.strptime(s, "%b").month
    else:
        try:
            return datetime.strptime(s, "%B").month
        except ValueError:
            return None


def lubridate(string, delim=',|-|/| '):
    """
    Guess the date format of a string, starting with ISO format.
    :param string:  str, serialization of a calendar date in unknown format
    :param delim:  str, a regular expression comprising all potential characters
                   delimiting date fields, separated by "|"
    :return:  datetime.date object
    """
    if string == "":
        return None

    try:
        # first try ISO format
        return date.fromisoformat(string)
    except ValueError:
        items = re.split(delim, string)

        if len(items) == 1:
            if string.isnumeric():
                year = str2year(string)
                if year is None:
                    return None
                return date(year, 7, 2)  # midpoint of year
            else:
                return None  # what is this even

        elif len(items) == 2:
            # try some month-year combos
            if items[0].isalpha() and items[1].isnumeric():
                return date(str2year(items[1]), str2month(items[0]), 15)
            elif items[0].isnumeric() and items[1].isalpha():
                return date(str2year(items[0]), str2month(items[1]), 15)
            elif items[0].isnumeric() and items[1].isnumeric():
                if len(items[0]) == 4 and len(items[1]) == 2:  # Ym
                    return date(int(items[0]), int(items[1]), 15)
                elif len(items[0]) == 2 and len(items[1]) == 4:  # mY
                    return date(int(items[1]), int(items[0]), 15)
                else:
                    return None  # ambiguous

        elif len(items) == 3:
            # not ISO format, first locate named month
            is_alpha = [x.isalpha() for x in items]
            if sum(is_alpha) != 1:
                return None
            month = str2month(items[is_alpha.index(True)])

            # next locate 4-digit year
            is_year = [x.isnumeric() and len(x) == 4 for x in items]
            if sum(is_year) != 1:
                return None
            year = str2year(items[is_year.index(True)])
            is_day = [not (is_alpha[i] or is_year[i]) for i in range(len(items))]
            day = int(items[is_day.index(True)])
            return date(year, month, day)

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
    parser.add_argument("--debug", action="store_true",
                        help="If set, write original and parsed dates to stdout "
                             "instead of relabeled FASTA.")
    parser.add_argument("--keep-all", action="store_true",
                        help="Write all sequences even if dates fail to parse.")
    parser.add_argument("--year", action="store_true",
                        help="Write only years of sample collection.")
    args = parser.parse_args()

    # load metadata and parse collection dates
    coldates = {}
    for row in DictReader(args.metadata):
        dt = lubridate(row[args.field])
        if args.debug:
            print(f"{row[args.field]},{dt}")
        coldates.update({row["header"]: dt})

    if args.debug:
        sys.exit()

    # apply parsed dates to sequence labels
    for record in SeqIO.parse(args.infile, args.format):
        header = record.name
        try:
            dt = coldates[header]
        except KeyError:
            print(f"ERROR: Failed to retrieve date for {header}")
            sys.exit()

        if dt is None and not args.keep_all:
            continue  # debug mode, quit without writing output

        if args.year:
            res = f">{record.description}_{dt.year}\n{record.seq}\n"
        else:
            res = f">{record.description}_{dt.isoformat()}\n{record.seq}\n"
        sys.stdout.write(res)
