"""Drop duplicate sequences from a given FASTA file.
"""
import argparse
import Bio.SeqIO


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", help="sequences to be deduplicated")
    parser.add_argument("--output", help="deduplicated sequences")

    args = parser.parse_args()

    sequences = Bio.SeqIO.parse(args.sequences, "fasta")
    observed_sequence_ids = set()

    with open(args.output, "w") as oh:
        for sequence in sequences:
            if sequence.id not in observed_sequence_ids:
                observed_sequence_ids.add(sequence.id)

                Bio.SeqIO.write(sequence, oh, "fasta")
            else:
                print(f"Skipping duplicate of {sequence.id}")
