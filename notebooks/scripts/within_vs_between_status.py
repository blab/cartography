import argparse
from augur.utils import read_node_data
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys

from Helpers import get_euclidean_data_frame

if __name__ == "__main__":

    # Initialize parsers
    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--embedding", required=True, help="the path to a dataframe csv file")
    parser.add_argument("--clades", required=True, help="a path to the clade status of the different strains in the build")
    parser.add_argument("--method", required=True, choices = ["pca", "mds", "t-sne", "umap"], help="the embedding used")
    parser.add_argument("--embedding-columns", nargs=2, required=True, help="list of the two columns to use as coordinates from the embedding data frame")
    parser.add_argument("--scaled-distance", required=True, help="path to distances between every strain")
    parser.add_argument("--output-figure", help="path for outputting as a PNG")
    parser.add_argument("--output-dataframe", help="path for outputting as a dataframe")
    
    args = parser.parse_args()
    
    #Error Handling
    if args.output_figure is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)

    #read in embedding file     
    embedding_df = pd.read_csv(args.embedding, index_col=0)
    #read in node data
    
    node_data = read_node_data(args.clades)
    
    sequences_by_name = list(embedding_df.index)

    sequences_list = []
    for sequence in sequences_by_name:
        if sequence in node_data["nodes"]:
            sequences_list.append(sequence)
        
    # Build a data frame of clade annotations per strain in the same order
    # as the sequences and the subsequent distance matrix.
    clade_annotations = pd.DataFrame([
        {"strain": sequence_name, "clade_membership": node_data["nodes"][sequence_name]["clade_membership"]}
        for sequence_name in sequences_list
    ])

    merged_df = embedding_df.merge(clade_annotations, on="strain")

    KDE_df = get_euclidean_data_frame(merged_df, args.embedding_columns[0], args.embedding_columns[1], args.method)


    KDE_df["scaled_distance"] = pd.read_csv(args.scaled_distance, index_col=0)

    if args.output_dataframe is not None:
        KDE_df.to_csv(args.output_dataframe)

    if args.output_figure is not None:
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax = sns.kdeplot(KDE_df.query("clade_status == 'within'")["scaled_distance"], label="Same clade", ax=ax)
        ax = sns.kdeplot(KDE_df.query("clade_status == 'between'")["scaled_distance"], label="Different clade", ax=ax)

        ax.set_xlabel("Scaled Euclidean distance from embedding")
        ax.set_ylabel("KDE density")

        fig.suptitle(args.method + ' KDE Plot', fontsize=16)
        sns.despine()

        plt.savefig(args.output_figure)