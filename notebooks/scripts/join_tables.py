"""Merge two tables.
"""
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", help="left table in join")
    parser.add_argument("--right", help="right table in join")
    parser.add_argument("--on", help="column to join on")
    parser.add_argument("--output", help="joined table")

    args = parser.parse_args()

    # Read tables.
    left = pd.read_csv(args.left, sep=None, engine="python")
    right = pd.read_csv(args.right, sep=None, engine="python")

    joined = left.merge(
        right,
        on=args.on,
    )

    joined.to_csv(
        args.output,
        sep="\t",
        index=False,
        header=True,
    )
