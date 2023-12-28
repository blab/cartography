#!/usr/bin/env python3
"""Calculate internal node embedding positions based on embeddings of observed
sequences and the topology of a given phylogenetic tree.
"""
import argparse
from collections import defaultdict

from augur.utils import read_tree
import numpy as np
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree", required=True, help="Newick file representing the phylogenetic relationship of samples in the given embedding with named internal nodes.")
    parser.add_argument("--embedding", required=True, help="CSV file with an embedding and cluster label per sample ('strain') for observed sequences.")
    parser.add_argument("--output", required=True, help="CSV file with the input embedding positions plus the calculated positions for the named internal nodes from the given tree.")

    args = parser.parse_args()

    # Load the tree.
    tree = read_tree(args.tree)

    # Load the embedding.
    embedding = pd.read_csv(args.embedding, index_col=0)

    # Find embedding position columns, ignoring all labels from clustering.
    embedding_columns = [
        column
        for column in embedding.columns
        if not "label" in column
    ]

    # Calculate positions for internal nodes.
    position_by_node = defaultdict(dict)
    for node in tree.find_clades(order="postorder"):
        for embedding_column in embedding_columns:
            if node.is_terminal():
                position_by_node[node.name][embedding_column] = embedding.at[
                    node.name,
                    embedding_column
                ]
            else:
                position_by_node[node.name][embedding_column] = np.mean([
                    position_by_node[child.name][embedding_column]
                    for child in node.clades
                ])

    # Create a data frame of positions per internal node.
    internal_node_positions = []
    for node in tree.find_clades(terminal=False, order="preorder"):
        internal_node = {
            embedding.index.name: node.name,
        }
        internal_node.update(position_by_node[node.name])
        internal_node_positions.append(internal_node)

    internal_node_positions = pd.DataFrame.from_records(
        internal_node_positions,
        index=embedding.index.name,
    )

    # Combine observed and internal node positions into a single data frame and
    # export it.
    pd.concat([embedding, internal_node_positions]).to_csv(
        args.output,
        header=True,
        index=True,
    )
