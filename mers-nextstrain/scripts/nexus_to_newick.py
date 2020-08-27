"""Convert a given NEXUS file to a Newick file.
"""
import argparse
import Bio.Phylo


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="a NEXUS tree file with sequences attached to the tree")
    parser.add_argument("output", help="a Newick file")

    args = parser.parse_args()

    tree = Bio.Phylo.read(args.tree, "nexus")

    for node in tree.find_clades():
        if node.name:
            node.name = node.name.strip("'").split("|")[0]

    Bio.Phylo.write(tree, args.output, "newick")
