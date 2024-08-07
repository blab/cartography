import argparse
import pandas as pd
import numpy as np


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", nargs="+", help="KDE plots")
    parser.add_argument("--output-csv", help="comma or tab-delimited file of attributes per node of the given tree")
    parser.add_argument("--output-table", help="rmarkdown table in file")
    parser.add_argument("--separator", default=',', help="separator between columns in the given tables")
    parser.add_argument("--disease-names", nargs="+", help="names of diseases, separated by space")
    parser.add_argument("--genetic-group-types", nargs="+", help="names of genetic group types per disease name listed above")

    args = parser.parse_args()

    columns_to_use = [
        "method",
        "normalized_vi",
        "n_predicted_clusters",
        "distance_threshold",
    ]
    embedding_name_by_abbreviation = {
        "pca": "PCA",
        "mds": "MDS",
        "t-sne": "t-SNE",
        "umap": "UMAP",
        "genetic": "genetic",
    }

    # Concatenate tables.
    tables = []
    for table_file, pathogen, genetic_group_type in zip(args.tables, args.disease_names, args.genetic_group_types):
        df = pd.read_csv(
            table_file,
            sep=args.separator,
        )

        if "distance_threshold" in df.columns:
            df = df.loc[:, columns_to_use].copy()
        else:
            df = df.loc[:, [column for column in columns_to_use if column != "distance_threshold"]].copy()
            df["distance_threshold"] = ""

        df["Pathogen Dataset"] = [pathogen] + [""] * (len(embedding_name_by_abbreviation) - 1)
        df["Genetic Group Type"] = [genetic_group_type] + [""] * (len(embedding_name_by_abbreviation) - 1)
        tables.append(df)

    df = pd.concat(tables, ignore_index=True)

    df["method"] = df["method"].map(embedding_name_by_abbreviation)
    df = df.rename(columns={
        "method": "Method",
        "normalized_vi": "Variation of Information (VI)",
        "n_predicted_clusters": "Number of clusters",
        "distance_threshold": "Threshold",
    })

    if args.output_csv:
        df.to_csv(args.output_csv, sep=args.separator, header=True)

    if args.output_table:
        latex_table = df.reset_index().to_latex(
            bold_rows=True,
            index=False,
            columns=[
                "Pathogen Dataset",
                "Genetic Group Type",
                "Method",
                "Number of clusters",
                "Variation of Information (VI)",
                "Threshold",
            ],
        )
        latex_table = latex_table.replace("bottomrule", "botrule")

        with open(args.output_table, "w", encoding="utf-8") as oh:
            oh.write(latex_table)
