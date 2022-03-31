import argparse
from augur.utils import read_node_data
import json
from numbers import Number


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Auspice JSON with a tree to annotate")
    parser.add_argument("--node-data", nargs="+", help="node data JSONs to annotate on the given tree.")
    parser.add_argument("--output", help="Auspice JSON annotated by the given node data")
    args = parser.parse_args()

    # Load the tree JSON.
    with open(args.tree, "r", encoding="utf-8") as fh:
        tree_json = json.load(fh)

    # Load the node data.
    node_data = read_node_data(args.node_data)
    node_data = node_data["nodes"]

    # Update the tree to include node data per node. Perform a preorder
    # traversal of the tree by processing each node first followed by its
    # immediate children. Perform this traversal iteratively by storing
    # references to original node objects in the JSON dictionary instead of
    # recursing through the tree and maintaining a large stack of function
    # calls.
    nodes_to_process = [tree_json["tree"]]
    type_by_key = {}
    while nodes_to_process:
        # Get the next node to process from the end of the stack.
        node = nodes_to_process.pop()
        print(f"Processing {node['name']}, {len(nodes_to_process)} nodes remain")

        # Annotate node, if data exist for it. While node data are key/value
        # pairs per strain, the Auspice JSON stores them as key/dict pairs where
        # the dict contains a key of "value" that indexes the value of the
        # key/value pair.
        if node["name"] in node_data:
            for key, value in node_data[node["name"]].items():
                node["node_attrs"][key] = {"value": value}

                if key not in type_by_key:
                    value_type = "categorical"
                    if isinstance(value, Number) and not "label" in key:
                        value_type = "continuous"

                    type_by_key[key] = value_type

        # Update list of nodes to process, if this node has children.
        if "children" in node:
            for child in node["children"]:
                nodes_to_process.append(child)

    # Update the metadata to include colorings for the given node data.
    for key, value_type in type_by_key.items():
        tree_json["meta"]["colorings"].append({
            "key": key,
            "title": key.replace("_", " "),
            "type": value_type,
        })

        # If the value type is categorical, add the key to the list of filters,
        # too.
        if value_type == "categorical":
            tree_json["meta"]["filters"].append(key)

    # Write out the annotated Auspice tree JSON.
    with open(args.output, "w", encoding="utf-8") as oh:
        json.dump(tree_json, oh)
