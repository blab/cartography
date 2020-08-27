"""Convert a given NEXUS file with sequences into a FASTA file.
"""
import argparse
import Bio.SeqIO


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="a NEXUS tree file with sequences attached to the tree")
    parser.add_argument("output", help="a FASTA file of the sequences attached to the tree")

    args = parser.parse_args()

    sequences = Bio.SeqIO.parse(args.tree, "nexus")

    Bio.SeqIO.write(sequences, args.output, "fasta")
