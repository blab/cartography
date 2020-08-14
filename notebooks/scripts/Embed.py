import argparse
from augur.utils import write_json
import Bio.SeqIO
from collections import OrderedDict
import pandas as pd
from scipy.spatial.distance import squareform
from sklearn.manifold import TSNE, MDS
from umap import UMAP
import sys

from Helpers import get_hamming_distances

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(dest="command")
    
    tsne = subparsers.add_parser("t-sne")
    tsne.add_argument("--alignment", required=True, help="an aligned FASTA file from a multiple sequence alignment")
    tsne.add_argument("--perplexity", default=30.0, type=float, help="the perplexity value for the tsne embedding")
    tsne.add_argument("--learning-rate", default=200.0, type=float, help="the learning rate value for the tsne embedding")
    tsne.add_argument("--output-node-data", help="outputting a node data JSON file")
    tsne.add_argument("--output-dataframe", help="outputting a csv file")
    
    umap = subparsers.add_parser("umap")
    umap.add_argument("--alignment", required=True, help="an aligned FASTA file from a multiple sequence alignment")
    umap.add_argument("--nearest-neighbors", default=200, type=int, help="the nearest neighbors value for the umap embedding")
    umap.add_argument("--min-dist", default=.5, type=float, help="the minimum distance value for the umap embedding")
    umap.add_argument("--output-node-data", help="outputting a node data JSON file")
    umap.add_argument("--output-dataframe", help="outputting a csv file")
    
    mds = subparsers.add_parser("mds")
    mds.add_argument("--alignment", required=True, help="an aligned FASTA file from a multiple sequence alignment")
    mds.add_argument("--components", default=10, type=int, help="the number of components for MDS")
    mds.add_argument("--output-node-data", help="outputting a node data JSON file")
    mds.add_argument("--output-dataframe", help="outputting a csv file")
    
    args = parser.parse_args()
    
    if args.output_node_data is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)
    options = ["pca", "mds", "tsne", "umap"]
    if args.command not in options:
        print("You must specify one of the options allowed:" + str(options[i] for i in options), file=sys.stderr)
        sys.exit(1)
    # Load alignment
    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(args.alignment, "fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)
    
    sequence_names = list(sequences_by_name.keys())
    
    # Calculate Distance Matrix
    hamming_distances = get_hamming_distances(
        sequences_by_name.values()
    )
    distance_matrix = squareform(hamming_distances)
    
    # Calculate Embedding
    if args.command == "t-sne":
        embedding_class = TSNE #conditionally learn name of subparser 
        embedding_parameters = {
            "metric": "precomputed",
            "perplexity": args.perplexity,
            "learning_rate": args.learning_rate
        }
    elif args.command == "umap":
        embedding_class = UMAP
        embedding_parameters = {
            "n_neighbors": args.nearest_neighbors,
            "min_dist": args.min_dist,
            "n_components": 2,
            "init": "spectral"
        }
    elif args.command == "mds":
        embedding_class = MDS
        embedding_parameters = {
            "dissimilarity": "precomputed",
            "n_components": args.components
        }

    embedder = embedding_class(**embedding_parameters)
    embedding = embedder.fit_transform(distance_matrix)
    
    print(embedding)
    
    # Output Embedding
        # create dictionary to be "wrapped" by write_json
        
    embedding_df = pd.DataFrame(embedding)
    if args.command == "mds" or args.command == "pca":
        embedding_df.columns=[args.command + str(i) for i in range(1,args.components + 1)]
    else:
        embedding_df.columns = [args.command.replace('-', '') + "_x" , args.command.replace('-', '') + "_y"]
        
    embedding_df.index = sequence_names
    if args.output_node_data is not None:
        embedding_dict = embedding_df.transpose().to_dict()
        write_json(embedding_dict,"results/embed_" + args.command.replace('-', '') + ".json")
    if args.output_dataframe is not None:
        embedding_df.to_csv("results/embed_" + args.command.replace('-', '') + ".csv")
    