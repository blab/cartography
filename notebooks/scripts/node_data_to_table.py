import sys
sys.path.append("../")

import argparse
from augur.utils import read_node_data, read_tree
import pandas as pd

from Helpers import get_y_positions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, help="Newick tree")
    parser.add_argument("--node-data", nargs="+", required=True, help="node data JSON(s) to extract attributes from")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data from internal nodes in output")
    parser.add_argument("--attributes", nargs="+", help="names of attributes to export from the given tree")
    parser.add_argument("--output", required=True, help="tab-delimited file of attributes per node of the given tree")

    args = parser.parse_args()

    # Load tree.
    tree = read_tree(args.tree)

    # Load node data.
    node_data = read_node_data(args.node_data)["nodes"]

    # Collect attributes per node from the tree to export.
    records = []

    heights = get_y_positions(tree)
    for node in tree.find_clades():
        if node.is_terminal() or args.include_internal_nodes:
            record = {
                "strain": node.name,
                "y_value": heights[node],
                "is_internal_node" : node.is_terminal()
            }

            for attribute in args.attributes:
                if attribute in node_data[node.name]:
                    record[attribute] = node_data[node.name][attribute]
                else:
                    print(f"Attribute '{attribute}' missing from node '{node.name}'", file=sys.stderr)

            records.append(record)

   # Convert records to a data frame and save as a tab-delimited file.
    df = pd.DataFrame(records)
    df.to_csv(args.output, sep="\t", header=True, index=False)
