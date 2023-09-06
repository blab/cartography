#!/usr/bin/env python3
import argparse
import pandas as pd
from scipy.spatial.distance import pdist, squareform


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--lineages", required=True, help="CSV table of parental lineages A and B and their recombinant offspring X")
    parser.add_argument("--method", required=True, help="name of the embedding method represented by the given embedding file. Use to annotate the output table.")
    parser.add_argument("--embedding", required=True, help="CSV of an embedding to calculate Euclidean distances from")
    parser.add_argument("--metadata", required=True, help="TSV of metadata containing strain name and lineage information required to find strains from recombinant lineages")
    parser.add_argument("--lineage-column", default="Nextclade_pango", help="name of the metadata column for the lineage information per strain")
    parser.add_argument("--output", required=True, help="CSV table of lineages with method and average pairwise distances between strains in each lineage")

    args = parser.parse_args()

    # Load lineages.
    lineages = pd.read_csv(args.lineages)

    # Load embedding.
    embedding = pd.read_csv(
        args.embedding,
        index_col="strain",
    )

    # Load metadata.
    metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        usecols=["strain", args.lineage_column],
        index_col="strain",
    )

    # Subset metadata to records present in the embedding.
    metadata_subset = metadata.loc[embedding.index]

    # Calculate average pairwise distance between lineages for each set of lineages.
    lineage_distances = []
    for row_id, (parental_A, parental_B, recombinant_X) in lineages.iterrows():
        # Find strains per lineage.
        parental_A_strains = set(metadata_subset.loc[metadata_subset[args.lineage_column] == parental_A].index.values)
        parental_B_strains = set(metadata_subset.loc[metadata_subset[args.lineage_column] == parental_B].index.values)
        recombinant_X_strains = set(metadata_subset.loc[metadata_subset[args.lineage_column] == recombinant_X].index.values)

        # Subset embedding to strains of interest.
        embedding_subset = embedding.loc[
            list(parental_A_strains | parental_B_strains | recombinant_X_strains),
        ]

        # Calculate all pairwise distances for the subset.
        distances = squareform(pdist(embedding_subset))

        # Calculate average pairwise distance between groups.
        A_index = [embedding_subset.index.get_loc(strain) for strain in parental_A_strains]
        B_index = [embedding_subset.index.get_loc(strain) for strain in parental_B_strains]
        X_index = [embedding_subset.index.get_loc(strain) for strain in recombinant_X_strains]
        A_B_distances = []
        A_X_distances = []
        B_X_distances = []

        for a in A_index:
            for b in B_index:
                A_B_distances.append(distances[a, b])

            for x in X_index:
                A_X_distances.append(distances[a, x])

        for b in B_index:
            for x in X_index:
                B_X_distances.append(distances[b, x])

        if len(A_B_distances) > 0 and len(A_X_distances) > 0 and len(B_X_distances) > 0:
            average_A_B_distance = sum(A_B_distances) / len(A_B_distances)
            average_A_X_distance = sum(A_X_distances) / len(A_X_distances)
            average_B_X_distance = sum(B_X_distances) / len(B_X_distances)

            lineage_distances.append({
                "parental_A": parental_A,
                "parental_B": parental_B,
                "recombinant_X": recombinant_X,
                "method": args.method,
                "distance_A_B": average_A_B_distance,
                "distance_A_X": average_A_X_distance,
                "distance_B_X": average_B_X_distance,
                "X_maps_closer_to_both_parentals": (
                    (average_A_X_distance < average_A_B_distance) and
                    (average_B_X_distance < average_A_B_distance)
                ),
                "X_maps_closer_to_any_parental": (
                    (average_A_X_distance < average_A_B_distance) or
                    (average_B_X_distance < average_A_B_distance)
                ),
            })
        else:
            print(f"Could not find enough samples to compare all three lineages: {parental_A}, {parental_B}, and {recombinant_X}")

    # Save average pairwise distances to a table.
    pd.DataFrame(lineage_distances).to_csv(
        args.output,
        index=False,
        float_format="%.4f",
    )
