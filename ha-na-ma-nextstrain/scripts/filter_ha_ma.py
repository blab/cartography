import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", nargs=3, help="sequences to filter, the one to filter first and one to filter by second")
    parser.add_argument("--output_fasta", nargs=3, help="FASTA files of split genomes")

    args = parser.parse_args()

    # Index sequences without loading them into memory. This gives us access to
    # names of strains in both files that we can cross-check.
    sequences_a = SeqIO.index(args.sequence[0], "fasta")
    sequences_b = SeqIO.index(args.sequence[1], "fasta")
    sequences_c = SeqIO.index(args.sequence[2], "fasta")

    # Identify shared sequences between the two sets.
    shared_strains = set(sequences_a.keys()) & set(sequences_b.keys() & set(sequences_c.keys()))
    print(f"Found {len(shared_strains)} between input sequence files.")

    # Write out shared strains for sequence set a.
    SeqIO.write(
        (sequences_a[strain] for strain in shared_strains),
        args.output_fasta[0],
        "fasta"
    )

    # Write out shared strains for sequence set b.
    SeqIO.write(
        (sequences_b[strain] for strain in shared_strains),
        args.output_fasta[1],
        "fasta"
    )

    SeqIO.write(
        (sequences_c[strain] for strain in shared_strains),
        args.output_fasta[2],
        "fasta"
    )
