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

    for node in tree.find_clades(terminal=True):
        if node.is_terminal() or args.include_internal_nodes:
            record = {
                "name": node.name
            }

            for attribute in attributes:
                if attribute in node.node_attrs:
                    record[attribute] = node.node_attrs[attribute]["value"]
                elif attribute in node.branch_attrs:
                    record[attribute] = node.branch_attrs[attribute]["value"]
                else:
                    print(f"Attribute '{attribute}' missing from node '{node.name}'", file=sys.stderr)

            records.append(record)

    tree_records = []
    heights = get_y_positions(tree)
    for node in tree.find_clades(terminal=True):
        tree_records.append(dict(strain=node.name, y=heights[node]))
    tree_records_df = pd.DataFrame(tree_records)
    # Convert records to a data frame and save as a tab-delimited file.
    df = pd.DataFrame(records)
    df.columns = ["strain"] + list(attributes)
    df["y_value"] = tree_records_df["y"]

    df.to_csv(args.output, sep="\t", header=True, index=False, float_format="%.2f")
