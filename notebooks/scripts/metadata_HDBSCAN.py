import argparse
import Bio.SeqIO
from collections import OrderedDict
import hdbscan
import numpy as np
import pandas as pd
import re
from sklearn.metrics import confusion_matrix, matthews_corrcoef
import sys

from Helpers import get_hamming_distances, get_euclidean_data_frame, get_embedding_columns_by_method


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--method", required=True, choices = ["pca", "mds", "t-sne", "umap"], help="the embedding used")
    parser.add_argument("--embedding", help="a csv file with embedding information for method")
    parser.add_argument("--metadata", help="clades information")
    parser.add_argument("--cluster-data", help="cluster data from embedding and assign labels given via HDBSCAN")
    parser.add_argument("--cluster-threshold", type=float, help="cluster data from embedding and assign labels given via HDBSCAN. Pass in a threshold.")
    parser.add_argument("--output", help="outputting a csv file of metadata info for HDBSCAN results")

    args = parser.parse_args()

    if args.output is not None:
        embedding_df = pd.read_csv(args.embedding)
        if args.metadata is not None:
            metadata_df = pd.read_csv(args.metadata, sep="\t")
            embedding_df = embedding_df.merge(metadata_df[["strain", "clade_membership"]], on="strain")
            KDE_df_cluster = get_euclidean_data_frame(sampled_df=embedding_df, column_for_analysis=f"{args.method}_label", embedding=args.method, column_list=get_embedding_columns_by_method(args.method))
            KDE_df_normal = get_euclidean_data_frame(sampled_df=embedding_df,column_for_analysis="clade_membership", embedding=args.method, column_list=get_embedding_columns_by_method(args.method))

            confusion_matrix_val = confusion_matrix(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
            matthews_cc_val = matthews_corrcoef(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
            if args.cluster_threshold is not None:
                output_df = pd.DataFrame([[args.method, matthews_cc_val, float(args.cluster_threshold), confusion_matrix_val[0][0], confusion_matrix_val[1][0], confusion_matrix_val[1][1], confusion_matrix_val[0][1]]], columns=["embedding", "MCC", "threshold", "TN", "FN", "TP", "FP"]).round(3)
                output_df.to_csv(args.output)
            if args.cluster_data is not None:
                max_df = pd.read_csv(args.cluster_data)
                output_df = pd.DataFrame([[args.method, matthews_cc_val, float(max_df['distance_threshold'].values.tolist()[0]), confusion_matrix_val[0][0], confusion_matrix_val[1][0], confusion_matrix_val[1][1], confusion_matrix_val[0][1]]], columns=["embedding", "MCC", "threshold", "TN", "FN", "TP", "FP"]).round(3)
                output_df.to_csv(args.output, index=False)
