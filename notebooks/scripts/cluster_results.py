import argparse
import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True, choices = ["pca", "mds", "t-sne", "umap"], help="the embedding used")
    parser.add_argument("--embedding", help="csv file with the embedding")
    parser.add_argument("--columns", nargs="+", help="the columns which the pdist will be calculated on")
    parser.add_argument("--threshold-information", help="the distance threshold csv to be used on HDBSCAN. if not provided, it will run without.")
    parser.add_argument("--column-threshold", help="name of threshold column in the threshold dataframe")
    parser.add_argument("--output", help="the csv path where the label of the data and strain name will be saved.")
    parser.add_argument("--output-figure", help="PNG with the results displayed graphically")

    args = parser.parse_args()

    embedding_df = pd.read_csv(args.embedding, index_col=0)
    
    if(args.threshold_information is not None):
        threshold_df = pd.read_csv(args.threshold_information)
        distance_threshold = threshold_df.loc[threshold_df['embedding'] == args.method][args.column_threshold].values.tolist()[0]
        sd = np.std(pdist(embedding_df))
        mean = np.mean(pdist(embedding_df))
        distance_threshold = (distance_threshold * sd) + mean
        distance_threshold = float(distance_threshold)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=15, cluster_selection_epsilon=distance_threshold)
        clusterer.fit(embedding_df)
        embedding_df[f"{args.method}_label"] = clusterer.labels_.astype(str)
        if(args.output):
            embedding_df.to_csv(args.output)
    else:
        clusterer = hdbscan.HDBSCAN(min_cluster_size=15)
        clusterer.fit(embedding_df)
        embedding_df[f"{args.method}_label"] = clusterer.labels_.astype(str)
        if(args.output):
            embedding_df.to_csv(args.output)