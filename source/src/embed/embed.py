import argparse

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


def get_hamming_distances(genomes):
    """Calculate pairwise Hamming distances between the given list of genomes
    and return the nonredundant array of values for use with scipy's squareform function.
    Bases other than standard nucleotides (A, T, C, G) are ignored.

    Parameters
    ----------
    genomes : list
        a list of strings corresponding to genomes that should be compared

    Returns
    -------
    list
        a list of distinct Hamming distances as a vector-form distance vector


    >>> genomes = ["ATGCT", "ATGCT", "ACGCT"]
    >>> get_hamming_distances(genomes)
    [0, 1, 1]
    >>> genomes = ["AT-GCT", "AT--CT", "AC--CT"]
    >>> get_hamming_distances(genomes)
    [0, 1, 1]

    """
    # Define an array of valid nucleotides to use in pairwise distance calculations.
    # Using a numpy array of byte strings allows us to apply numpy.isin later.

    nucleotides = np.array([b'A', b'T', b'C', b'G', b'a', b't', b'c', b'g']) # take into account lowercase letters as well as uppercase

    # Convert genome strings into numpy arrays to enable vectorized comparisons.
    genome_arrays = [
        np.frombuffer(genome.encode(), dtype="S1")
        for genome in genomes
    ]

    # Precalculate positions of valid bases (A, T, C, and G) in each genome to speed up later comparisons.
    valid_bases = [
        np.isin(genome_array, nucleotides)
        for genome_array in genome_arrays
    ]

    # Calculate Hamming distance between all distinct pairs of genomes at valid bases.
    # The resulting list is a reduced representation of a symmetric matrix that can be
    # converted to a square matrix with scipy's squareform function:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
    hamming_distances = []
    for i in range(len(genomes)):
        # Only compare the current genome, i, with all later genomes.
        # This avoids repeating comparisons or comparing each genome to itself.
        for j in range(i + 1, len(genomes)):
            # Find all mismatches between these two genomes.
            mismatches = genome_arrays[i] != genome_arrays[j]

            # Count the number of mismatches where both genomes have valid bases.
            hamming_distances.append((mismatches & valid_bases[i] & valid_bases[j]).sum())

    return hamming_distances


def embed(args):
    # TODO: Create a default cluster distance threshold if none given.

    # Setting Random seed for numpy
    np.random.seed(seed=args.random_seed)

    if args.output_dataframe is None:
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
            "n_init": 2
        }

    if args.command != "pca":
        embedder = embedding_class(**embedding_parameters)
        embedding = embedder.fit_transform(distance_matrix)

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
