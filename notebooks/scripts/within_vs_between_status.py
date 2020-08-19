import argparse
from augur.utils import read_node_data
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

if __name__ == "__main__":

    # Initialize parsers
    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--embedding-csv", required=True, help="the path to a dataframe csv file")
    parser.add_argument("--node-data", required=True, help="a path to the clade status of the different strains in the build")
    parser.add_argument("--method", required=True, help="the embedding used")
    parser.add_argument("--output-PNG", help="path for outputting as a PNG")
    parser.add_argument("--output-dataframe", help="path for outputting as a dataframe")
    
    args = parser.parse_args()
    
    #Error Handling
    if args.output_PNG is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)
        
    options = ["pca", "mds", "t-sne", "umap"] #add genetic here? 
    if args.method not in options:
        print("You must specify one of the options allowed:" + str(options[i] for i in options), file=sys.stderr)
        sys.exit(1)
    #read in embedding file     
    embedding_df = pd.read_csv(args.embedding_csv, index_col=0)
    #read in node data
    
    node_data = read_node_data(clades_path)
    
    sequences_by_name = embedding_df.index

    for sequence in sequences_by_name:
        if sequence not in node_data["nodes"]:
            sequences_by_name = sequences_by_name.remove(sequence)
        
    # Build a data frame of clade annotations per strain in the same order
    # as the sequences and the subsequent distance matrix.
    clade_annotations = pd.DataFrame([
        {"strain": sequence_name, "clade": node_data["nodes"][sequence_name]["clade_membership"]}
        for sequence_name in sequence_names
    ])
    
    def assign_clade_status_to_pairs(clade_annotations):
        """Assign clade status to all pairs in the given list of indices and the given data frame of clade annotations.
        
        Outputs a vector in condensed distance matrix format such that all nonredundant pairs of strains are represented.
        
        """
        clade_statuses = []
        for i in range(len(clade_annotations)):
            for j in range(i + 1, len(clade_annotations)):
                same_clade = clade_annotations.loc[i, "clade"] == clade_annotations.loc[j, "clade"]
                clade_statuses.append(int(same_clade))
                
        return np.array(clade_statuses)
    
    
    clade_annotations["clade_status"] = assign_clade_status_to_pairs(clade_annotations)

    if args.output_PNG is not None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax = sns.kdeplot(clade_annotations.query("clade_status == 1")["clade_status"], label="Same clade", ax=ax)
        ax = sns.kdeplot(clade_annotations.query("clade_status == 0")["clade_status"], label="Different clade", ax=ax)

        ax.set_xlabel("Scaled Euclidean distance from embedding")
        ax.set_ylabel("KDE density")

        fig.suptitle(args.method + 'KDE Plot', fontsize=16)
        sns.despine()