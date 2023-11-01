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
            # Deduplicate by the sequence description, since many strain names
            # contain spaces and do not get parsed properly. We know the
            # description contains the full defline of each FASTA record and
            # that each line can contain pipe-delimited metadata with the strain
            # name in the first field.
            sequence_id = sequence.description.split("|")[0]

            if sequence_id not in observed_sequence_ids:
                observed_sequence_ids.add(sequence_id)

                Bio.SeqIO.write(sequence, oh, "fasta")
            else:
                print(f"Skipping duplicate of {sequence_id}")
