import sys
sys.path.append("../")

import numpy as np
import argparse
from augur.io import read_metadata
from augur.utils import read_node_data, read_tree, annotate_parents_for_tree
import pandas as pd

from Helpers import get_y_positions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, help="Newick tree")
    parser.add_argument("--metadata", help="metadata with attributes to include in output table. When a given attribute exists in metadata and node data, this script will select the node data value over the metadata value.")
    parser.add_argument("--node-data", nargs="+", required=True, help="node data JSON(s) to extract attributes from")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data from internal nodes in output")
    parser.add_argument("--attributes", nargs="+", help="names of attributes to export from the given tree")
    parser.add_argument("--mutation-length-attribute", default="mutation_length", help="attribute name of the mutation length value stored in the branch lengths node data")
    parser.add_argument("--output", required=True, help="tab-delimited file of attributes per node of the given tree")

    args = parser.parse_args()

    # Load tree.
    tree = read_tree(args.tree)
    tree = annotate_parents_for_tree(tree)

    # Load node data.
    node_data = read_node_data(args.node_data)["nodes"]

    # Load metadata.
    metadata = {}
    if args.metadata:
        metadata = read_metadata(args.metadata)
        metadata = metadata.to_dict(orient="index")

    if "divergence" in args.attributes:
        # loop through mutation lengths, making them cumulative and store in "divergence"
        # look through branch lengths data for the mutation length attribute (e.g., "mutation_length")
        # start with root of tree
        for node in tree.find_clades():
            if getattr(node, "parent") is None:
                node_data[node.name]["divergence"] = node_data[node.name][args.mutation_length_attribute]
            else:
                node_data[node.name]["divergence"] = node_data[node.parent.name]["divergence"] + node_data[node.name][args.mutation_length_attribute]

    # Collect attributes per node from the tree to export.
    records = []

    heights = get_y_positions(tree)
    for node in tree.find_clades():
        if node.is_terminal() or args.include_internal_nodes:
            record = {
                "strain": node.name,
                "y_value": heights[node],
                "parent_name": node.parent.name if getattr(node, "parent") else "",
                "is_internal_node" : not node.is_terminal()
            }
            


            # Try to load attribute values from node data first and then metadata.
            for attribute in args.attributes:
                if attribute in node_data[node.name]:
                    record[attribute] = node_data[node.name][attribute]
                elif attribute in metadata.get(node.name, {}):
                    record[attribute] = metadata[node.name][attribute]
                else:
                    print(f"Attribute '{attribute}' missing from node '{node.name}'", file=sys.stderr)

            records.append(record)

   # Convert records to a data frame and save as a tab-delimited file.
    df = pd.DataFrame(records)

    def get_parent_y(row):
        parent_name = row['parent_name']
        if parent_name is not None:
            parent_y = df.loc[df['strain'] == parent_name, 'y_value'].values
            if len(parent_y) > 0:
                return parent_y[0]
            else:
                return row['y_value']

        return np.nan

    def get_parent_mutation_length(row):
        parent_name = row['parent_name']
        if parent_name is not None:
            parent_y = df.loc[df['strain'] == parent_name, 'divergence'].values
            if len(parent_y) > 0:
                return parent_y[0]
            else:
                return row['divergence']

        return np.nan

    if "divergence" in args.attributes:
        df['parent_y'] = df.apply(get_parent_y, axis=1)
        df['parent_mutation'] = df.apply(get_parent_mutation_length, axis=1)

    df.to_csv(args.output, sep="\t", header=True, index=False)
