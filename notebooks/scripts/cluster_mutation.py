from augur.io import read_sequences
from augur.utils import read_node_data
from collections import OrderedDict
import pandas as pd
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--reference-sequence", required=True, help="FASTA file of the reference sequence used for the alignment")
    parser.add_argument("--alignment", required=True, help="aligned FASTA file of diseases")
    parser.add_argument("--metadata", required=True, help="metadata with clade information")
    parser.add_argument("--metadata-column", default="MCC", help="metadata column to find clade information")
    parser.add_argument("--valid-characters", nargs="+", default=["A", "C", "T", "G", "-"], help="list of valid characters to consider in pairwise comparisons with the reference")
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

    print(f"Ignoring leading and trailing gaps in the reference at sites: {ignored_sites}")
    print(f"Valid characters: {args.valid_characters}")
    sequences_by_name = OrderedDict()

    for sequence in read_sequences(args.alignment):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = list(sequences_by_name.keys())

    # Index cluster name by sequence name (Python dict) from metadata
    # NOTE TO JOHN: Not sure if json is necessary here since it's the default MCC value given
    if args.metadata.endswith(".json"):
        node_data = read_node_data(args.metadata)
        clade_annotations = pd.DataFrame([
            {"strain": strain, "clade": annotations[args.metadata_column]}
            for strain, annotations in node_data["nodes"].items()
            if strain in sequences_by_name
        ])
        clade_annotations = clade_annotations.set_index("strain")
        clade_annotations = clade_annotations[clade_annotations.clade != "unassigned"]
    elif args.metadata.endswith(".csv") :
        clade_annotations = pd.read_csv(args.metadata)
        clade_annotations["clade"] = clade_annotations[args.metadata_column]
        clade_annotations = clade_annotations[clade_annotations.clade != -1]
        clade_annotations = clade_annotations[["strain", "clade"]]
        clade_annotations = clade_annotations.set_index("strain")

    # Find mutations per cluster relative to reference as Python dictionary of sets indexed by cluster name
    clade = clade_annotations.groupby(["clade"])
    strains_per_cluster_reference = {}
    for clade in clade_annotations.groupby(["clade"]):
        strains_per_cluster = set(clade[1].index)
        strains_per_cluster_reference[clade[0]] = strains_per_cluster

    mutations_per_cluster_reference = {}
    for clade in strains_per_cluster_reference:
        mutations=[]
        for name in strains_per_cluster_reference[clade]:
            strain = sequences_by_name[name]
            # compare each strain to the reference
            for site in range(len(reference)):
                if all((
                    site not in ignored_sites,
                    strain[site] != reference[site],
                    strain[site] in args.valid_characters,
                    reference[site] in args.valid_characters,
                )):
                    # Report the site in 1-based coordinates.
                    mutations.append(
                        {
                            "cluster": name,
                            "site": site + 1,
                            "reference_allele": reference[site],
                            "alternate_allele": strain[site],
                        }
                    )

        mutations = pd.DataFrame(mutations)
        mutations["mutation"] = mutations.apply(
            lambda row: str(row["site"]) + row["alternate_allele"],
            axis=1,
        )
        mutations_per_cluster_reference[clade] = set(mutations["mutation"])

    # Remove mutations that are present in all clusters (reference-only mutations)
    # Intersection of all sets in sets by cluster name to find shared mutations
    shared_reference_mutations = mutations_per_cluster_reference[list(mutations_per_cluster_reference)[0]]
    for mutation in mutations_per_cluster_reference:
        shared_reference_mutations = shared_reference_mutations.intersection(mutations_per_cluster_reference[mutation])

    # Remove resulting intersection from each cluster set
    for mutation in mutations_per_cluster_reference:
        mutations_per_cluster_reference[mutation] = mutations_per_cluster_reference[mutation] - shared_reference_mutations

    # For each pair of clusters, find mutations that only exist in “from” and “to” clusters
    # Nested for loop of clusters, excluding self-comparisons
    mutations = []
    for mutation_1 in mutations_per_cluster_reference:
        for mutation_2 in mutations_per_cluster_reference:
            if (mutation_1 != mutation_2):
                # Find set difference for from and to (from - to)
                set_diff = mutations_per_cluster_reference[mutation_1] - mutations_per_cluster_reference[mutation_2]
                # Report from cluster, to cluster, mutation (e.g., 750G)
                for item in set_diff:
                    mutations.append(
                        { "from": mutation_1,
                        "to": mutation_2,
                        "mutation": item
                        }
                    )
    mutations_df = pd.DataFrame(mutations)
    mutations_df.to_csv(args.output, index=False)
