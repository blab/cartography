"""Root a given tree, prune the root, and collapse low support internal nodes.
"""
import argparse
import Bio.Phylo


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="sequences to be deduplicated")
    parser.add_argument("--root", help="name of sample to root tree with and then remove")
    parser.add_argument("--minimum-node-support", default=0, type=int, help="minimum bootstrap support for an internal node below which the node will be collapsed")
    parser.add_argument("--output", help="rooted tree with low-support nodes collapsed")

    args = parser.parse_args()

    tree = Bio.Phylo.read(args.tree, "newick")

    # Root tree with the requested sample.
    tree.root_with_outgroup(args.root)

    # Prune the root from the tree.
    tree.prune(args.root)

    # Collapse nodes with low support.
    print(f"Started with {len(tree.get_nonterminals())} internal nodes")
    tree.collapse_all(
        lambda c: c.confidence is not None and c.confidence < args.minimum_node_support
    )
    print(f"Ended with {len(tree.get_nonterminals())} internal nodes")

    # Ladderize.
    tree.ladderize()

    with open(args.output, "w", encoding="utf-8") as oh:
        Bio.Phylo.write(tree, oh, "newick")
