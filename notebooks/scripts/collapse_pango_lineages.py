#!/usr/bin/env python3
import argparse
from augur.io import read_metadata
from collections import defaultdict
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
    count_by_lineage = pango_counts.to_dict(defaultdict(int))

    # Sort lineages by count in ascending order so we collapse the smaller
    # lineages into their parents first, allowing larger small lineages to get
    # populated first so they don't need to be collapsed.
    low_count_lineages = pango_counts[
        pango_counts < args.min_samples
    ].sort_values().to_dict()

    new_lineage_by_old = {}
    for lineage, count in low_count_lineages.items():
        # Check whether the current lineage is still below the threshold, given
        # that earlier iterations could have incremented its counts above and we
        # no longer want to collapse this lineage.
        if count_by_lineage.get(lineage, 0) >= args.min_samples:
            print(f"{lineage} is now big enough to skip collapsing (N={count_by_lineage[lineage]})")
            continue

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
            # If there is no parent, then we either try to collapse the current
            # lineage into the highest level parent we could reach or we keep
            # the lineage as it is.
            if lineage in new_lineage_by_old:
                print(f"Ran out of parents for {lineage}, adding its {count_by_lineage[lineage]} samples to {new_lineage_by_old[lineage]}")
                parent = new_lineage_by_old[lineage]
                count_by_lineage[parent] += count_by_lineage.pop(lineage)
            else:
                print(f"Could not find any parents for {lineage}")

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
