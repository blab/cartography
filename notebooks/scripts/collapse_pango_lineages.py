#!/usr/bin/env python3
import argparse
from augur.io import read_metadata
from pango_aliasor.aliasor import Aliasor


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--metadata", required=True, help="metadata TSV with pango lineages to be collapsed")
    parser.add_argument("--pango-column", default="Nextclade_pango", help="name of the column in the metadata with pango lineages to collapse")
    parser.add_argument("--min-samples", type=int, default=10, help="minimum number of samples required for a lineage to not be collapsed")
    parser.add_argument("--new-column", default="Nextclade_pango_collapsed", help="name of the new column with collapsed pango lineages")
    parser.add_argument("--output", required=True, help="metadata TSV with new column representing collapsed pango lineages")

    args = parser.parse_args()

    aliasor = Aliasor()
    metadata = read_metadata(args.metadata)

    pango_counts = metadata[args.pango_column].value_counts()
    count_by_lineage = pango_counts.to_dict()
    low_count_lineages = pango_counts[pango_counts < args.min_samples].to_dict()

    new_lineage_by_old = {}
    for lineage, count in low_count_lineages.items():
        parent = aliasor.parent(lineage)

        # When we reach the top of the pango hierarchy, the parent will be an
        # empty string.
        while parent:
            if count_by_lineage.get(parent, 0) + count >= args.min_samples:
                # If adding the current lineage counts to its parent makes this
                # a high count lineage, map the lineage to its parent, remove it
                # from counts, and update parent counts.
                print(f"Adding {count} samples from {lineage} to its parent {parent}")
                new_lineage_by_old[lineage] = parent
                del count_by_lineage[lineage]
                count_by_lineage[parent] += count
                break
            else:
                # If we can't add this lineage to an existing parent entry, keep
                # the mapping of lineage to its parent and try the next parent
                # up the hierarchy.
                print(f"Moving from {lineage} to {parent}")
                new_lineage_by_old[lineage] = parent
                parent = aliasor.parent(parent)

        if parent == "":
            # If there is no parent, move on to the next lineage, keeping counts
            # for this lineage as they are.
            print(f"Ran out of parents for {lineage}")
            new_lineage_by_old[lineage] = lineage

    # Map lineages, keeping the old if new not defined.
    metadata[args.new_column] = metadata[args.pango_column].apply(
        lambda lineage: new_lineage_by_old.get(lineage, lineage)
    )

    # Save new metadata.
    metadata.to_csv(
        args.output,
        sep="\t",
        index=True,
    )