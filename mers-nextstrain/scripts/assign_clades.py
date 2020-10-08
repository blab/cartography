#!/usr/bin/env python
# coding: utf-8
import argparse
from augur.utils import read_metadata, read_tree, write_json


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", help="tab-delimited metadata")
    parser.add_argument("--tree", help="Newick tree with internal node names")
    parser.add_argument("--output", help="node data JSON with clade membership annotations")

    args = parser.parse_args()

    metadata, metadata_fields = read_metadata(args.metadata)
    tree = read_tree(args.tree)

    # Look for clades for which all children have the same host. To do this,
    # make a postorder traversal of the tree such that each internal node gets
    # marked with the host of its children if all children have the same host.
    # Otherwise, the internal node is marked with a host of `None` to note that
    # its children were sampled from multiple hosts.
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            node.host = metadata[node.name]["host"]
        else:
            # Find all unique hosts of this node's children.
            hosts = list({child.host for child in node.clades})
            if len(hosts) == 1 and hosts[0] is not None:
                node.host = hosts[0]
            else:
                node.host = None

    # Create a "node data" JSON file for auspice with an annotation for "clade
    # membership". Assign clade names using an incremental count
    clades = {}
    clade_count = 0

    for node in tree.find_clades(terminal=False, order="preorder"):
        # Only assign clade names to internal nodes with a defined host.
        if node.host is not None and node.name not in clades:
            clades[node.name] = {
                "clade_membership": f"clade_{clade_count}"
            }

            # All descendants of this node should belong to the same host,
            # if our host assignment algorithm above worked properly.
            for child in node.find_clades():
                clades[child.name] = {
                    "clade_membership": f"clade_{clade_count}"
                }

            # The next internal node will belong to a different clade.
            clade_count += 1

    # Write out the node data JSON for consumption by `augur export`.
    write_json({"nodes": clades}, args.output)
