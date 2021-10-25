import argparse
from augur.utils import write_json
import Bio.SeqIO
from collections import OrderedDict
import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from scipy.spatial.distance import squareform, pdist
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS
import sys
from umap import UMAP

from Helpers import get_hamming_distances, get_euclidean_data_frame


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--distance-matrix", help="a csv distance matrix that can be read in by pandas, index column as row 0")
    parser.add_argument("--alignment", help="an aligned FASTA file to create a distance matrix with")
    parser.add_argument("--cluster-data", help="cluster data from embedding and assign labels given via HDBSCAN")
    parser.add_argument("--cluster-threshold", type=float, help="cluster data from embedding and assign labels given via HDBSCAN. Pass in a threshold.")
    parser.add_argument("--random-seed", default = 314159, type=int, help="an integer used as the random seed for reproducible results")
    parser.add_argument("--output-node-data", help="outputting a node data JSON file")
    parser.add_argument("--output-dataframe", help="outputting a csv file")
    parser.add_argument("--output-figure", help="plot of the embedding, for debugging purposes")

    #parser.add_argument("--method-params" help="csv file from grid search")
        #if method params exists and command line exists, command line overrides it
    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    pca = subparsers.add_parser("pca")
    pca.add_argument("--components", default=10, type=int, help="the number of components for PCA")
    pca.add_argument("--explained-variance", default="results/explained_variance_pca.png", help="the path for the explained variance table")

    tsne = subparsers.add_parser("t-sne")
    tsne.add_argument("--perplexity", default=30.0, type=float, help="the perplexity value for the tsne embedding")
    tsne.add_argument("--learning-rate", default=200.0, type=float, help="the learning rate value for the tsne embedding")

    umap = subparsers.add_parser("umap")
    umap.add_argument("--nearest-neighbors", default=200, type=int, help="the nearest neighbors value for the umap embedding")
    umap.add_argument("--min-dist", default=.5, type=float, help="the minimum distance value for the umap embedding")

    mds = subparsers.add_parser("mds")
    mds.add_argument("--components", default=10, type=int, help="the number of components for MDS")

    args = parser.parse_args()
    # Checking that the input fits the restrictions

    # Setting Random seed for numpy
    np.random.seed(seed=args.random_seed)

    if args.output_node_data is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)

    if args.alignment is None and args.command == "pca":
        print("You must specify an alignment for pca, not a distance matrix", file=sys.stderr)
        sys.exit(1)

    # getting or creating the distance matrix

    if args.distance_matrix is not None:
        distance_matrix  = pd.read_csv(args.distance_matrix, index_col=0)

    elif args.alignment is not None:
        sequences_by_name = OrderedDict()

        for sequence in Bio.SeqIO.parse(args.alignment, "fasta"):
            sequences_by_name[sequence.id] = str(sequence.seq)

        sequence_names = list(sequences_by_name.keys())
        if args.command != "pca":
            # Calculate Distance Matrix

            hamming_distances = get_hamming_distances(
                sequences_by_name.values()
            )
            distance_matrix = pd.DataFrame(squareform(hamming_distances))
            distance_matrix.index = sequence_names

    # Calculate Embedding

    if args.command == "pca":
        sequences_by_name = OrderedDict()

        for sequence in Bio.SeqIO.parse(args.alignment, "fasta"):
            sequences_by_name[sequence.id] = str(sequence.seq)

        sequence_names = list(sequences_by_name.keys())

        numbers = list(sequences_by_name.values())[:]
        for i in range(0,len(list(sequences_by_name.values()))):
            numbers[i] = re.sub(r'[^AGCT]', '5', numbers[i])
            numbers[i] = list(numbers[i].replace('A','1').replace('G','2').replace('C', '3').replace('T','4'))
            numbers[i] = [int(j) for j in numbers[i]]

        genomes_df = pd.DataFrame(numbers)
        genomes_df.columns = ["Site " + str(k) for k in range(0,len(numbers[i]))]


        #performing PCA on my pandas dataframe
        pca = PCA(n_components=args.components,svd_solver='full') #can specify n, since with no prior knowledge, I use None
        principalComponents = pca.fit_transform(genomes_df)

        # Create a data frame from the PCA embedding.
        embedding = principalComponents
        embedding_df = pd.DataFrame(principalComponents)
        embedding_df.index = sequence_names

    if args.command == "t-sne":
        embedding_class = TSNE
        embedding_parameters = {
            "metric": "precomputed",
            "perplexity": args.perplexity,
            "learning_rate": args.learning_rate,
            "random_state" : args.random_seed,
            "square_distances": True,
        }
    elif args.command == "umap":
        embedding_class = UMAP
        embedding_parameters = {
            "n_neighbors": args.nearest_neighbors,
            "min_dist": args.min_dist,
            "n_components": 2,
            "init": "spectral",
            "random_state" : args.random_seed
        }
    elif args.command == "mds":
        embedding_class = MDS
        embedding_parameters = {
            "dissimilarity": "precomputed",
            "n_components": args.components,
            "n_jobs": 1,
            "n_init": 2,
            "random_state": args.random_seed
        }

    if args.command != "pca":
        embedder = embedding_class(**embedding_parameters)
        embedding = embedder.fit_transform(distance_matrix)

        print(embedding)

        # Output Embedding
            # create dictionary to be "wrapped" by write_json

        embedding_df = pd.DataFrame(embedding)
        embedding_df.index = list(distance_matrix.index)

    if args.command == "mds" or args.command == "pca":
        embedding_df.columns=[args.command + str(i) for i in range(1,args.components + 1)]
    else:
        embedding_df.columns = [args.command.replace('-', '') + "_x" , args.command.replace('-', '') + "_y"]

    if args.command == "pca":

        #add explained variance as the first row of the dataframe
        explained_variance = pd.DataFrame([round(pca.explained_variance_ratio_[i],4) for i in range(0,len(pca.explained_variance_ratio_))], columns=["explained variance"])
        explained_variance["principal components"] = [i for i in range(1, args.components + 1)]
        explained_variance.to_csv(args.explained_variance, index=False)

    clusterer = None
    if args.cluster_threshold is not None:
        cluster = float(args.cluster_threshold)
        clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=float(cluster))
    elif args.cluster_data is not None:
        max_df = pd.read_csv(args.cluster_data)
        clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=float(max_df.where(max_df["method"] == args.command).dropna(subset = ['distance_threshold'])[["distance_threshold"]].values.tolist()[0][0]))

    if clusterer is not None:
        clusterer_default = hdbscan.HDBSCAN()
        clusterer.fit(embedding_df)
        clusterer_default.fit(embedding_df)
        embedding_df[f"{args.command}_label"] = clusterer.labels_.astype(str)
        embedding_df[f"{args.command}_label_default"] = clusterer_default.labels_.astype(str)

    if args.output_node_data is not None:
        embedding_dict = embedding_df.transpose().to_dict()
        write_json({"nodes": embedding_dict}, args.output_node_data)

    if args.output_dataframe is not None:
        embedding_df.to_csv(args.output_dataframe, index_label="strain")

    if args.output_figure:
        plot_data = {
            "x": embedding[:, 0],
            "y": embedding[:, 1],
        }

        if clusterer is not None:
            plot_data["cluster"] = clusterer.labels_.astype(str)
        else:
            plot_data["cluster"] = "0"

        plot_df = pd.DataFrame(plot_data)
        ax = sns.scatterplot(
            data=plot_df,
            x="x",
            y="y",
            hue="cluster",
            alpha=0.5,
        )
        plt.savefig(args.output_figure)
        plt.close()
