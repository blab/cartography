import argparse
from augur.utils import write_json
import pandas as pd
import sys


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--metadata", help="a decompressed tsv metadata file that can be read into pandas")
    parser.add_argument("--output", help="a clades.json file to be used by the KDE plots")
    parser.add_argument("--col-name", help="cluster data from embedding and assign labels given via HDBSCAN")
    
    args = parser.parse_args()

    metadata_df = pd.read_csv(args.metadata, sep="\t", index_col=0)

    metadata_df.rename(columns={args.col_name:"clade_membership"}, inplace=True)
    clades_df = metadata_df[["clade_membership"]]
    
    if args.output is not None:
        clades_dict = clades_df.transpose().to_dict()
        write_json({"nodes": clades_dict}, args.output)