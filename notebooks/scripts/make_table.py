"""Merge two or more tables as data frames.
"""
import argparse
from functools import reduce
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", nargs="+", help="tables to concatenate")
    parser.add_argument("--separator", default="\t", help="separator between columns in the given tables")
    parser.add_argument("--suffixes", nargs=2, help="what to add when two columns have the same value")
    parser.add_argument("--output", help="concatenated table")

    args = parser.parse_args()

    # Read tables.
    tables = []
    for i in range(0, len(args.tables)):
        tables.append(pd.read_csv(args.tables[i], sep=args.separator, dtype=str))

    if args.suffixes is not None:
        df = reduce(lambda x, y: pd.merge(x, y, on = 'strain', suffixes=(args.suffixes[0], args.suffixes[1])), tables)
    else:
        df = reduce(lambda x, y: pd.merge(x, y, on = 'strain'), tables)


    df.to_csv(args.output, sep=args.separator, header=True, index=False)
