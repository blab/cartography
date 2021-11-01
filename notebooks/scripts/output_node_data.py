"""Create node json (augur) from pandas dataframe.
"""
import argparse
import pandas as pd
from augur.utils import write_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", help="table to make node data from")
    parser.add_argument("--separator", default=",", help="separator between columns in the given tables")
    parser.add_argument("--node_name", default="nodes", help="what to name the node value in the auspice json")
    parser.add_argument("--output", help="json file")

    args = parser.parse_args()

    if args.output is not None:
        embedding_dict = pd.read_csv(
            args.table,
            sep=args.separator,
            index_col=0
        ).transpose().to_dict()
        write_json({args.node_name: embedding_dict}, args.output)
