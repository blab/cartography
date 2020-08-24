import argparse
from augur.utils import write_json
import numpy as np
import pandas as pd
import re
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS
import sys
from umap import UMAP

from Helpers import get_hamming_distances


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--distance-matrix", required=True, help="a csv distance matrix that can be read in by pandas, index column as row 0")
    parser.add_argument("--output-node-data", help="outputting a node data JSON file")
    parser.add_argument("--output-dataframe", help="outputting a csv file")

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    pca = subparsers.add_parser("pca")
    pca.add_argument("--components", default=10, type=int, help="the number of components for PCA")

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

    if args.output_node_data is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)

    distance_matrix  = pd.read_csv(args.distance_matrix, index_col=0)

    # Calculate Embedding
    if args.command == "pca":
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
        embedding_df = pd.DataFrame(principalComponents)
        
    if args.command == "t-sne":
        embedding_class = TSNE  
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

    if args.command != "pca":  
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
        
    embedding_df.index = list(distance_matrix.index)
    if args.command == "pca":
    
        #add explained variance as the first row of the dataframe
        explained_variance = pd.DataFrame([round(pca.explained_variance_ratio_[i],4) for i in range(0,len(pca.explained_variance_ratio_))], columns=["explained variance"])
        explained_variance["principal components"] = [i for i in range(0, args.components)] 
        explained_variance.to_csv("results/explained_variance_pca.csv", index=False)
        
    if args.output_node_data is not None:
        embedding_dict = embedding_df.transpose().to_dict()
        write_json(embedding_dict,args.output_node_data)
    if args.output_dataframe is not None:
        embedding_df.to_csv(args.output_dataframe)
    