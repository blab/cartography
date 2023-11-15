#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from functools import reduce
import operator


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--alignment", required=True, help="alignment from which to drop all columns with N characters")
    parser.add_argument("--output", required=True, help="alignment without N characters")
    args = parser.parse_args()

    # Load alignment.
    alignment = AlignIO.read(args.alignment, "fasta")

    # Find positions to drop from the alignment.
    drop_indexes = [
        i
        for i in range(alignment.get_alignment_length())
        if "N" in alignment[:, i] or "-" in alignment[:, i]
    ]
    print(f"Dropping {len(drop_indexes)} from alignment of length {alignment.get_alignment_length()}.")

    # Collect positions that will not be dropped.
    edited_alignment = reduce(
        operator.add,
        (
            alignment[:, index:index+1]
            for index in range(alignment.get_alignment_length())
            if index not in drop_indexes
        )
    )

    # Save the new alignment.
    AlignIO.write(
        edited_alignment,
        args.output,
        "fasta",
    )
