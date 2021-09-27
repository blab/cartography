"""
pathogen-embed.

Reduced dimension embeddings for pathogen sequences.
"""

__version__ = "0.0.2"
__author__ = 'Sravani Nanduri, John Huddleston'
__credits__ = 'Bedford Lab, Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA'

import argparse
import sys 
from .embed import embed

def make_parser():
    parser = argparse.ArgumentParser(description = "Reduced dimension embeddings for pathogen sequences", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--distance-matrix", help="a distance matrix that can be read in by pandas, index column as row 0")
    parser.add_argument("--separator", default=",", help="separator between columns in the given distance matrix")
    parser.add_argument("--alignment", help="an aligned FASTA file to create a distance matrix with. Make sure the strain order in this file matches the order in the distance matrix.")
    parser.add_argument("--cluster-data", help="The file (same separator as distance-matrix) that contains the distance threshold by which to cluster data in the embedding and assign labels given via HDBSCAN (https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html). If no value is given in cluster-data or cluster-threshold, the default distance threshold of 0.0 will be used.")
    parser.add_argument("--cluster-threshold", type=float, help="The float value for the distance threshold by which to cluster data in the embedding and assign labels via HDBSCAN. If no value is given in cluster-data or cluster-threshold, the default distance threshold of 0.0 will be used.")
    parser.add_argument("--random-seed", default = 314159, type=int, help="an integer used for reproducible results.")
    parser.add_argument("--output-dataframe", help="a csv file outputting the embedding with the strain name and its components.")
    parser.add_argument("--output-figure", help="outputs a PNG plot of the embedding")

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    pca = subparsers.add_parser("pca", description="Principal Component Analysis")
    pca.add_argument("--components", default=10, type=int, help="the number of components for PCA")
    pca.add_argument("--explained-variance", help="the path for the CSV explained variance for each component")

    tsne = subparsers.add_parser("t-sne", description="t-distributed Stochastic Neighborhood Embedding")
    tsne.add_argument("--perplexity", default=30.0, type=float, help="The perplexity is related to the number of nearest neighbors. Because of this, the size of the dataset is proportional to the best perplexity value (large dataset -> large perplexity). Values between 5 and 50 work best. The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")
    tsne.add_argument("--learning-rate", default=200.0, type=float, help="The learning rate for t-SNE is usually between 10.0 and 1000.0. Values out of these bounds may create innacurate results. The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")

    umap = subparsers.add_parser("umap", description="Uniform Manifold Approximation and Projection")
    umap.add_argument("--nearest-neighbors", default=200, type=int, help="Nearest neighbors controls how UMAP balances local versus global structure in the data (finer detail patterns versus global structure). This value is proportional to the size of the data (large dataset -> large nearest neighbors. The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")
    umap.add_argument("--min-dist", default=.5, type=float, help="Minimum Distance controls how tightly packed the UMAP embedding is. While it does not change the structure of the data, it does change the embedding's shape. The default value is the value consistently the best for pathogen analyses, results from an exhaustive grid search.")

    mds = subparsers.add_parser("mds", description="Multidimensional Scaling")
    mds.add_argument("--components", default=10, type=int, help="the number of components for MDS")

    return parser


def run(argv):
    args = make_parser().parse_args(argv)
    try:
        return embed(args)
    except Exception as error:
        print(error, file=sys.stderr)
        sys.exit(1)