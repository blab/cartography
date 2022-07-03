from augur.io import read_sequences
import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="the name of the aligned consensus file WITH reference")
    parser.add_argument("--output", required=True, help="outputting a tsv file of mutations per cluster")

    args = parser.parse_args()


    sequences = read_sequences(args.alignment)

    # The reference is the first sequence.
    reference = next(sequences)

    # The rest of the sequences are from the clusters.
    clusters = list(sequences)

    ignored_clusters = ["unassigned"]
    mutations = []
    for cluster in clusters:
        if cluster.name in ignored_clusters:
            continue

        for site in range(len(reference.seq)):
            if cluster.seq[site] != reference.seq[site]:
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