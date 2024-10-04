import argparse
import altair as alt
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True, help="path to TSV of counts by lineage and cluster")
    parser.add_argument("--lineage-column", required=True, help="name of column in counts file for lineages")
    parser.add_argument("--cluster-column", required=True, help="name of column in counts file for clusters")
    parser.add_argument("--output", required=True, help="PNG dot plot of counts by lineage and cluster")

    args = parser.parse_args()

    df = pd.read_csv(
        args.counts,
        sep="\t",
    )

    chart = alt.Chart(df).mark_circle().encode(
        x=alt.X(f"{args.lineage_column}:N", title="Pango lineage"),
        y=alt.Y(f"{args.cluster_column}:N", title="Cluster from t-SNE"),
        size="count:Q",
    ).properties(
        width=400,
        height=400,
    )
    chart.save(args.output, format="png", ppi=300)
