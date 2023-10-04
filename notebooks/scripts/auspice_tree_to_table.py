import sys
sys.path.append("../")

import argparse
from augur.utils import json_to_tree
import Bio
import Bio.Phylo
import json
import pandas as pd
import sys

from Helpers import get_y_positions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice tree JSON")
    parser.add_argument("output", help="tab-delimited file of attributes per node of the given tree")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data from internal nodes in output")
    parser.add_argument("--attributes", nargs="+", help="names of attributes to export from the given tree")

    args = parser.parse_args()

    # Load tree from JSON.
    with open(args.tree, "r") as fh:
        tree_json = json.load(fh)

    tree = json_to_tree(tree_json)

    # Collect attributes per node from the tree to export.
    records = []

    if args.attributes:
        attributes = args.attributes
    else:
        attributes = sorted(list(tree.root.node_attr.keys()) + list(tree.root.branch_attrs.keys()))

    heights = get_y_positions(tree)
    for node in tree.find_clades(terminal=True):
        if node.is_terminal() or args.include_internal_nodes:
            record = {
                "strain": node.name,
                "y_value": heights[node],
            }

            for attribute in attributes:
                if attribute in node.node_attrs:
                    # Most node attributes have a dictionary with their value
                    # stored by a "value" key, but some core attributes like
                    # "div" or "accession" are scalar values.
                    if type(node.node_attrs[attribute]) is dict and "value" in node.node_attrs[attribute]:
                        record[attribute] = node.node_attrs[attribute]["value"]
                    else:
                        record[attribute] = node.node_attrs[attribute]
                elif attribute in node.branch_attrs:
                    record[attribute] = node.branch_attrs[attribute]["value"]
                else:
                    print(f"Attribute '{attribute}' missing from node '{node.name}'", file=sys.stderr)

            records.append(record)

   # Convert records to a data frame and save as a tab-delimited file.
    df = pd.DataFrame(records)
    df.to_csv(args.output, sep="\t", header=True, index=False)
