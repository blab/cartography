"""
calculate the average and std dev genetic distance within and between groups given a distance matrix,
a per-strain "group" definition file (e.g., a TSV with Nextstrain clade or a CSV with HDBSCAN cluster
labels from a t-SNE embedding), and the column of the given file to use for groups.
"""
import sys
sys.path.append("../")

import argparse
import json
import numpy as np
import pandas as pd
from scipy import stats
import sys

from Helpers import get_euclidean_data_frame


def describe_dict(result):

    description_dict = {
        "mean": np.mean(result),
        "median" : np.median(result),
        "std" : np.std(result),
        "group" : args.group_column
    }
    return description_dict

if __name__ == "__main__":

    # Initialize parsers
    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # group column instead of cluster (naming)
    parser.add_argument("--distance-matrix", required=True, help="the path to a CSV genetic distance matrix")
    parser.add_argument("--metadata", required=True, help="a CSV/TSV path to a metadata table with the group status of the different strains in the build")
    parser.add_argument("--group-column", required=True, help="the column corresponding to grouping status")
    parser.add_argument("--ignored-group-labels", nargs="+", default=["-1", "unassigned"], help="group labels from group-column to ignore")
    parser.add_argument("--output", required=True,  help="path for outputting as a dataframe")

    args = parser.parse_args()

    distance_matrix = pd.read_csv(args.distance_matrix, index_col=0)

    node_data = pd.read_csv(args.metadata, sep="\t" if args.metadata.endswith(".tsv") else ",", index_col="strain")
    node_dict = node_data.transpose().to_dict()

    group_annotations = pd.DataFrame([
        {"strain": sequence, args.group_column: node_dict[sequence][args.group_column]}
        for sequence in node_data.index
    ])

    # ignoring -1 cluster grouping
    group_annotations = group_annotations[~group_annotations[args.group_column].isin(args.ignored_group_labels)]

    # difference between clades
    distance_matrix.columns = distance_matrix.index
    indices_to_drop = distance_matrix[~distance_matrix.index.isin(group_annotations["strain"])].dropna(how = 'all')
    distance_matrix = distance_matrix[distance_matrix.index.isin(group_annotations["strain"])].dropna(how = 'all')
    distance_matrix = distance_matrix.drop(indices_to_drop.index, axis=1)
    distance_matrix["strain"] = distance_matrix.index

    # difference between clusters
    merged_group_df = distance_matrix.merge(group_annotations, on="strain")
    distance_group_df = get_euclidean_data_frame(sampled_df=merged_group_df, column_for_analysis=args.group_column, embedding="genetic")

    # between cluster
    group_btw = describe_dict(np.array(distance_group_df[distance_group_df["clade_status"] == 0]["distance"]))

    # within cluster
    group_within = describe_dict(np.array(distance_group_df[distance_group_df["clade_status"] == 1]["distance"]))

    final_dict = {
        "between": group_btw,
        "within": group_within
    }

    output_df = pd.DataFrame.from_dict(final_dict, orient='index')
    output_df.index = output_df.index.rename("comparison")
    output_df.to_csv(args.output)
