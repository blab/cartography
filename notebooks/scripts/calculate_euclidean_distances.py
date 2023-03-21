#!/usr/bin/env python3
import argparse
import pandas as pd
from scipy.spatial.distance import pdist, squareform


if __name__ == '__main__':
    # Configure command line interface.
    parser = argparse.ArgumentParser()
    parser.add_argument("--embedding", required=True, help="embedding data frame")
    parser.add_argument("--output", required=True, help="CSV file of pairwise distances for each pair of tips in the embedding.")
    args = parser.parse_args()

    # Load embedding.
    embedding = pd.read_csv(
        args.embedding,
        index_col="strain",
    )

    # Drop label columns.
    label_columns = [column for column in embedding.columns if "label" in column]
    embedding = embedding.drop(columns=label_columns)

    # Calculate Euclidean distance.
    distances = pdist(embedding)

    # Save distances.
    tip_names = embedding.index.values

    distance_matrix = pd.DataFrame(squareform(distances), index=tip_names)
    pd.DataFrame(distance_matrix).to_csv(
        args.output,
        index=True,
        float_format="%.4f",
    )
