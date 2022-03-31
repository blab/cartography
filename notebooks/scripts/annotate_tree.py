import argparse
from augur.utils import write_json
from collections import OrderedDict
import pandas as pd
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Auspice JSON with a tree to annotate")
    parser.add_argument("--node-data", help="node data JSONs to annotate on the given tree.")
    parser.add_argument("--output", help="a clades.json file to be used by the KDE plots")
    parser.add_argument("--col-name", help="cluster data from embedding and assign labels given via HDBSCAN")
    args = parser.parse_args()

    if args.sequences is not None:
        sequences_by_name = OrderedDict()

        for sequence in Bio.SeqIO.parse(args.sequences, "fasta"):
            sequences_by_name[sequence.id] = str(sequence.seq)

        sequence_names_val = list(sequences_by_name.keys())
        print(len(sequence_names_val))
    metadata_df = pd.read_csv(args.metadata, sep="\t", index_col=0)

    if args.sequences is not None:
        metadata_df = metadata_df.loc[sequence_names_val]
        print(metadata_df)

    metadata_df.rename(columns={args.col_name:"clade_membership"}, inplace=True)
    clades_df = metadata_df[["clade_membership"]]

    if args.output is not None:
        clades_dict = clades_df.transpose().to_dict()
        write_json({"nodes": clades_dict}, args.output)
