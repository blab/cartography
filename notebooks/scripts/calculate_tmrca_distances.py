#!/usr/bin/env python3
import argparse
import Bio.Phylo
import json
import numpy as np
from scipy.spatial.distance import squareform


def annotate_parents_for_tree(tree):
    """Annotate each node in the given tree with its parent.
    >>> import io
    >>> tree = Bio.Phylo.read(io.StringIO("(A, (B, C))"), "newick")
    >>> not any([hasattr(node, "parent") for node in tree.find_clades()])
    True
    >>> tree = annotate_parents_for_tree(tree)
    >>> tree.root.parent is None
    True
    >>> all([hasattr(node, "parent") for node in tree.find_clades()])
    True
    """
    tree.root.parent = None
    for node in tree.find_clades(order="level"):
        for child in node.clades:
            child.parent = node

    # Return the tree.
    return tree


def calculate_distance_between_tips(tree, tip_a, tip_b, dates):
    """Calculate distance between tips.
    """
    # Find the MRCA of the given tips in the given tree.
    mrca = tip_a
    while tip_b.name not in mrca.terminals:
        mrca = mrca.parent

    # Calculate the distance from each tip to the MRCA.
    mrca_date = dates[mrca.name]["num_date"]
    tip_a_date = dates[tip_a.name]["num_date"]
    tip_b_date = dates[tip_b.name]["num_date"]
    latest_tip_date = max(tip_a_date, tip_b_date)

    # Return the maximum distance of the two tips.
    return latest_tip_date - mrca_date


if __name__ == '__main__':
    # Configure command line interface.
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, help="Newick tree")
    parser.add_argument("--dates", required=True, help="Augur node data JSON file with 'num_date' numeric date annotations per node (as from augur refine output)")
    parser.add_argument("--output", required=True, help="TSV file of pairwise distances to the tMRCA for each pair of tips in the tree, sorted by tip name")
    args = parser.parse_args()

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load dates.
    with open(args.dates, "r", encoding="utf-8") as fh:
        dates = json.load(fh)["nodes"]

    # Annotate tips per node to speed up MRCA calculations. Makes a single pass
    # through the tree in postorder to store a set of all terminals descending
    # from each node. This uses more memory, but it allows faster identification
    # of MRCAs between any pair of tips in the tree and speeds up pairwise
    # distance calculations by orders of magnitude.
    for node in tree.find_clades(order="postorder"):
        node.terminals = set()
        for child in node.clades:
            if child.is_terminal():
                node.terminals.add(child.name)
            else:
                node.terminals.update(child.terminals)

    # Calculate distances between all tips in condensed format.
    tips = sorted(tree.get_terminals(), key=lambda tip: tip.name)
    n_tips = len(tips)
    distances = np.zeros(n_tips * (n_tips - 1) // 2)
    for j in range(n_tips):
        for i in range(0, j):
            index = n_tips * i + j - ((i + 2) * (i + 1)) // 2
            distance = calculate_distance_between_tips(
                tree,
                tips[i],
                tips[j],
                dates,
            )
            distances[index] = distance

    # Save distances.
    np.savetxt(args.output, distances)
