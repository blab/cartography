import argparse
import numpy as np
import pandas as pd

from Helpers import variation_of_information


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--method", required=True, help="name of the embedding method")
    parser.add_argument("--true-clusters", required=True, help="metadata TSV or CSV with true cluster labels per strain")
    parser.add_argument("--true-clusters-column", required=True, help="column with true cluster labels in the given input table")
    parser.add_argument("--predicted-clusters", required=True, help="embedding CSV file with predicted cluster labels per strain")
    parser.add_argument("--predicted-clusters-column", required=True, help="column with predicted cluster labels in the given input table")
    parser.add_argument("--ignored-clusters", nargs="+", default=["-1", "unassigned"], help="list of cluster labels to ignore when calculating cluster accuracy")
    parser.add_argument("--output", required=True, help="CSV file with accuracy of clusters calculated by variation of information (VI) score")

    args = parser.parse_args()

    # Load true clusters.
    true_clusters = pd.read_csv(
        args.true_clusters,
        sep="\t" if args.true_clusters.endswith(".tsv") else ",",
        index_col="strain",
        usecols=["strain", args.true_clusters_column],
        dtype=str,
    )

    # Load predicted clusters.
    predicted_clusters = pd.read_csv(
        args.predicted_clusters,
        index_col="strain",
        usecols=["strain", args.predicted_clusters_column],
        dtype=str,
    )

    # Drop ignored cluster labels from both inputs.
    n_ignored_true_clusters = true_clusters[args.true_clusters_column].isin(args.ignored_clusters).sum()
    true_clusters = true_clusters[
        ~true_clusters[args.true_clusters_column].isin(args.ignored_clusters)
    ].copy()

    n_ignored_predicted_clusters = predicted_clusters[args.predicted_clusters_column].isin(args.ignored_clusters).sum()
    predicted_clusters = predicted_clusters[
        ~predicted_clusters[args.predicted_clusters_column].isin(args.ignored_clusters)
    ].copy()

    # Count the number of clusters assigned by HDBSCAN after dropping ignored
    # clusters. We also want to ignore the "-1" cluster label which represents
    # the lack of clustering.
    n_predicted_clusters = len(set(predicted_clusters[args.predicted_clusters_column].drop_duplicates().values) - {"-1"})

    # Join true and predicted labels by strain name, keeping only records that
    # appear in both sets.
    clusters = true_clusters.join(predicted_clusters, how="inner").reset_index()

    # Get sets of strains per true and predicted cluster.
    true_values = [
        list(cluster["strain"].values)
        for cluster_label, cluster in clusters.groupby(args.true_clusters_column)
    ]

    predicted_values = [
        list(cluster["strain"].values)
        for cluster_label, cluster in clusters.groupby(args.predicted_clusters_column)
    ]

    # Calculate variation of information (VI) score. Lower values indicate less
    # distance between true and predicted clusters.
    normalized_vi = variation_of_information(
        true_values,
        predicted_values,
        normalized=True,
    )

    output_df = pd.DataFrame([{
        "method": args.method,
        "predicted_clusters_column": args.predicted_clusters_column,
        "normalized_vi": np.round(normalized_vi, 2),
        "n_predicted_clusters": n_predicted_clusters,
        "n_predicted_cluster_samples": predicted_clusters.shape[0],
        "n_ignored_predicted_clusters": n_ignored_predicted_clusters,
        "n_true_cluster_samples": true_clusters.shape[0],
        "n_ignored_true_clusters": n_ignored_true_clusters,
        "n_vi_cluster_samples": clusters.shape[0],
    }])
    output_df.to_csv(args.output, index=False)
