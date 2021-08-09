#take in a json file, return a json file with the selected column renamed. 
import argparse
from augur.utils import read_node_data, write_json
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--clades", help="sequences to filter, the one to filter first and one to filter by second")
    parser.add_argument("--embedding", help="embedding dataframe from which the index can be used for the strains")
    parser.add_argument("--differentiator-column", help="the column name to change")
    parser.add_argument("--rename-column", help="name to rename the value to")
    parser.add_argument("--output", help="FASTA files of split genomes")

    args = parser.parse_args()
    embedding_df = pd.read_csv(args.embedding, index_col=0)
    #read in node data
    if args.clades is not None:
        node_data = read_node_data(args.clades)

        sequences_by_name = list(embedding_df.index)

        sequences_list = []
        for sequence in sequences_by_name:
            if sequence in node_data["nodes"]:
                sequences_list.append(sequence)

        # Build a data frame of clade annotations per strain in the same order
        # as the sequences and the subsequent distance matrix.
        clade_annotations = pd.DataFrame([
            {"strain": sequence_name, args.differentiator_column: node_data["nodes"][sequence_name][args.differentiator_column]}
            for sequence_name in sequences_list
        ])
        print(clade_annotations)
        clade_annotations = clade_annotations.rename(columns={args.differentiator_column: args.rename_column})
        clade_annotations.index = clade_annotations["strain"]
        clade_annotations = clade_annotations.drop(['strain'], axis=1)
        print(clade_annotations)
        if args.output is not None:
            embedding_dict = clade_annotations.transpose().to_dict()
            write_json({"nodes": embedding_dict}, args.output)