from augur.io import read_sequences
import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="the name of the aligned consensus file WITH reference")
    parser.add_argument("--ignored-characters", nargs="+", default=["X"], help="list of characters to ignore in pairwise comparisons with the reference")
    parser.add_argument("--output", required=True, help="outputting a tsv file of mutations per cluster")

    args = parser.parse_args()


    sequences = read_sequences(args.alignment)

    # The reference is the first sequence.
    reference = next(sequences)

    # Create a list of gap sites at the beginning and end of the reference
    # sequence to ignore from the alignment. This happens when the reference is
    # missing sequences that are present in the consensus sequences. These
    # differences usually reflect missing data in the original sequencing of the
    # reference and not a biologically-relevant insertion.
    ignored_sites = []
    site = 0
    while reference.seq[site] == "-" and site < len(reference.seq):
        ignored_sites.append(site)
        site += 1

    site = len(reference.seq) - 1
    while reference.seq[site] == "-" and site >= 0:
        ignored_sites.append(site)
        site -= 1

    print(f"Ignoring leading and trailing gaps in the reference at sites: {ignored_sites}")

    # The rest of the sequences are from the clusters.
    clusters = list(sequences)

    ignored_clusters = ["unassigned"]
    mutations = []
    for cluster in clusters:
        if cluster.name in ignored_clusters:
            continue

        for site in range(len(reference.seq)):
            # Find mismatches of characters we don't ignore at sites we don't
            # ignore.
            if (site not in ignored_sites and
                cluster.seq[site] != reference.seq[site] and
                cluster.seq[site] not in args.ignored_characters and
                reference.seq[site] not in args.ignored_characters):
                # Report the site in 1-based coordinates.
                mutations.append(
                    {
                        "cluster": cluster.name,
                        "site": site + 1,
                        "reference_allele": reference.seq[site],
                        "alternate_allele": cluster.seq[site],
                    }
                )

    mutations = pd.DataFrame(mutations)
    mutations["mutation"] = mutations.apply(
        lambda row: row["reference_allele"] + str(row["site"]) + row["alternate_allele"],
        axis=1,
    )

    mutations["is_cluster_specific"] = False
    for cluster_name, cluster in mutations.groupby("cluster"):
        mutations_in_cluster = set(cluster["mutation"].values)
        mutations_in_other_clusters = set(
            mutations.loc[
                mutations["cluster"] != cluster_name,
                "mutation"
            ].values
        )
        mutations_specific_to_cluster = mutations_in_cluster - mutations_in_other_clusters
        
        cluster_specific_index = (
            (mutations["cluster"] == cluster_name) & 
            mutations["mutation"].isin(mutations_specific_to_cluster)
        )
        mutations.loc[
            cluster_specific_index,
            "is_cluster_specific"
        ] = True

        print(f"{cluster_name}: {mutations_specific_to_cluster}")

    mutations.to_csv(
        args.output,
        sep="\t",
        index=False,
    )