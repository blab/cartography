import argparse
from augur.utils import write_json
import Bio.SeqIO
from collections import OrderedDict
from scipy.spatial.distance import squareform
from sklearn.manifold import TSNE
import sys

from Helpers import get_hamming_distances

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers()
    
    tsne = subparsers.add_parser("t-sne")
    tsne.add_argument("--alignment", required=True, help="an aligned FASTA file from a multiple sequence alignment")
    tsne.add_argument("--perplexity", default=30.0, type=float, help="the perplexity value for the tsne embedding")
    tsne.add_argument("--learning-rate", default=200.0, type=float, help="the learning rate value for the tsne embedding")
    tsne.add_argument("--output-node-data", help="outputting a node data JSON file")
    tsne.add_argument("--output-dataframe", help="outputting a csv file")
    
    args = parser.parse_args()
    
    if args.output_node_data is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
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
    embedding_class = TSNE #conditionally learn name of subparser 
    embedding_parameters = {
        "metric": "precomputed",
        "perplexity": args.perplexity,
        "learning_rate": args.learning_rate
    }
    embedder = embedding_class(**embedding_parameters)
    embedding = embedder.fit_transform(distance_matrix)
    
    print(embedding)
    # Output Embedding
        # create dictionary to be "wrapped" by write_json