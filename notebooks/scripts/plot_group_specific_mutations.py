#!/usr/bin/env python3
import argparse
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns


def rename_clade_name_field(cluster_column):
    name_by_method = {
        "pca": "PCA",
        "mds": "MDS",
        "t-sne": "t-SNE",
        "umap": "UMAP",
    }

    if cluster_column == "clade_membership":
        return "Nextstrain clade"
    else:
        new_name = cluster_column.replace("_label", "").replace("_", " ").replace(" for Nextstrain clade", "")
        return name_by_method.get(new_name, new_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mutation-table", required=True, help="table of mutations per pathogen dataset and genetic group")
    parser.add_argument("--output", required=True, help="plot of group-specific mutations per pathogen dataset and genetic group")

    args = parser.parse_args()

    df = pd.read_csv(args.mutation_table)

    filters = (
        (df["cluster_count"] == 1) &
        (df["dataset_name"] != "seasonal-flu-h3n2-ha-na-2016-2018") &
        (~df["cluster_column"].str.contains("pango"))
    )

    count_df = df[filters].groupby([
        "dataset_name",
        "cluster_column",
        "distinct_clusters",
    ]).aggregate(
        mutation_count=("mutation", "count"),
    ).reset_index()

    count_df["cluster_name"] = count_df["cluster_column"].map(rename_clade_name_field)

    g = sns.catplot(
        data=count_df,
        x="cluster_name",
        y="mutation_count",
        row="dataset_name",
        aspect=3.0,
        height=2,
        sharey=False,
        order=[
            "Nextstrain clade",
            "PCA",
            "MDS",
            "t-SNE",
            "UMAP",
            "genetic",
        ],
        row_order=[
            "seasonal-flu-h3n2-ha-2016-2018",
            "seasonal-flu-h3n2-ha-2018-2020",
            "sars-cov-2-2020-2022",
            "sars-cov-2-2022-2023",
        ],
        alpha=0.5,
    )

    g.set_axis_labels(
        x_var="genetic group",
        y_var="number of\ngroup-specific\nmutations",
    )
    for ax in g.axes.flatten():
        ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
