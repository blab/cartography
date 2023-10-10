import argparse
from augur.io import read_sequences
from augur.utils import read_node_data
from collections import OrderedDict
import pandas as pd
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--reference-sequence", required=True, help="FASTA file of the reference sequence used for the alignment")
    parser.add_argument("--alignment", required=True, help="aligned FASTA file of diseases")
    parser.add_argument("--metadata", required=True, help="metadata with clade information")
    parser.add_argument("--metadata-column", required=True, help="metadata column to find clade information")
    parser.add_argument("--ignored-clusters", nargs="+", default=["-1", "unassigned"], help="list of cluster labels to ignore when calculating cluster-specific mutations")
    parser.add_argument("--valid-characters", nargs="+", default=["A", "C", "T", "G", "-"], help="list of valid characters to consider in pairwise comparisons with the reference")
    parser.add_argument("--min-allele-count", type=int, default=0, help="minimum number of strains in a cluster with a given alternate allele required to include the allele in cluster-specific mutations")
    parser.add_argument("--min-allele-frequency", type=float, default=0.0, help="minimum frequency of an allele in a cluster to include allele in cluster-specific mutations")
    parser.add_argument("--output", help="the name of the csv file to be outputted")

    args = parser.parse_args()

    reference = str(next(read_sequences(args.reference_sequence)).seq)

    # Create a list of gap sites at the beginning and end of the reference
    # sequence to ignore from the alignment. This happens when the reference is
    # missing sequences that are present in the consensus sequences. These
    # differences usually reflect missing data in the original sequencing of the
    # reference and not a biologically-relevant insertion.
    ignored_sites = []
    site = 0
    while reference[site] == "-" and site < len(reference):
        ignored_sites.append(site)
        site += 1

    site = len(reference) - 1
    while reference[site] == "-" and site >= 0:
        ignored_sites.append(site)
        site -= 1

    print(f"Ignoring leading and trailing gaps in the reference at sites: {ignored_sites}", file=sys.stderr)
    print(f"Valid characters: {args.valid_characters}", file=sys.stderr)
    sequences_by_name = OrderedDict()

    for sequence in read_sequences(args.alignment):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = list(sequences_by_name.keys())

    # Index cluster name by sequence name from metadata.
    if args.metadata.endswith(".json"):
        node_data = read_node_data(args.metadata)
        clade_annotations = pd.DataFrame([
            {"strain": strain, "cluster": annotations[args.metadata_column]}
            for strain, annotations in node_data["nodes"].items()
            if strain in sequences_by_name
        ])
    elif args.metadata.endswith(".csv") :
        clade_annotations = pd.read_csv(args.metadata)
        clade_annotations["cluster"] = clade_annotations[args.metadata_column]
        clade_annotations = clade_annotations[["strain", "cluster"]]

    clade_annotations = clade_annotations.set_index("strain")

    # Remove records for cluster labels we ignore. For example, we often ignore
    # cluster labels that represent unclustered samples ("-1").
    clade_annotations["cluster"] = clade_annotations["cluster"].astype(str)
    clade_annotations = clade_annotations[
        ~clade_annotations["cluster"].isin(args.ignored_clusters)
    ].copy()

    # Find mutations per cluster relative to reference as Python dictionary of sets indexed by cluster name
    strains_per_cluster_reference = {}
    for clade in clade_annotations.groupby("cluster"):
        strains_per_cluster = set(clade[1].index)
        strains_per_cluster_reference[clade[0]] = strains_per_cluster

    all_mutation_counts = []
    for cluster, cluster_strains in strains_per_cluster_reference.items():
        mutations = []
        for strain in cluster_strains:
            strain_sequence = sequences_by_name[strain]

            # Compare each strain to the reference.
            for site in range(len(reference)):
                if all((
                    site not in ignored_sites,
                    strain_sequence[site] != reference[site],
                    strain_sequence[site] in args.valid_characters,
                    reference[site] in args.valid_characters,
                )):
                    # Report the site in 1-based coordinates.
                    mutations.append(
                        {
                            "cluster": cluster,
                            "strain": strain,
                            "site": site + 1,
                            "reference_allele": reference[site],
                            "alternate_allele": strain_sequence[site],
                        }
                    )

        mutations = pd.DataFrame(mutations)
        mutations["mutation"] = mutations.apply(
            lambda row: str(row["site"]) + row["alternate_allele"],
            axis=1,
        )

        mutation_counts = mutations.groupby("mutation")["strain"].count().reset_index().rename(columns={"strain": "count"})
        mutation_counts["frequency"] = mutation_counts["count"] / len(cluster_strains)

        # Filter by minimum allele count and frequency.
        mutation_counts = mutation_counts[mutation_counts["count"] >= args.min_allele_count].copy()
        mutation_counts = mutation_counts[mutation_counts["frequency"] >= args.min_allele_frequency].copy()

        mutation_counts["cluster"] = cluster

        all_mutation_counts.append(mutation_counts)

    # Count the number of clusters with each distinct alternate allele and find
    # the specific clusters with each allele.
    all_mutation_counts = pd.concat(all_mutation_counts, ignore_index=True)

    mutation_cluster_counts = all_mutation_counts.groupby(
        "mutation"
    ).aggregate(
        cluster_count=("cluster", "count"),
        distinct_clusters=("cluster", "unique")
    ).reset_index()

    # Remove mutations that are present in all clusters (reference-only mutations).
    total_clusters = all_mutation_counts["cluster"].drop_duplicates().shape[0]
    mutation_cluster_counts = mutation_cluster_counts[mutation_cluster_counts["cluster_count"] < total_clusters].copy()

    # Format distinct clusters as a comma-delimited list of names.
    mutation_cluster_counts["distinct_clusters"] = mutation_cluster_counts["distinct_clusters"].apply(
        lambda clusters: ",".join(clusters)
    )

    # Annotate column used for clusters.
    mutation_cluster_counts["metadata_column"] = args.metadata_column

    # Save mutations and their cluster information.
    mutation_cluster_counts.to_csv(args.output, index=False)
