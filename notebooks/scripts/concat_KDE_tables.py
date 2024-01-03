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

    # Concatenate tables.
    tables = []
    for table_file, pathogen, genetic_group_type in zip(args.tables, args.disease_names, args.genetic_group_types):
        df = pd.read_csv(table_file, sep=args.separator)
        df["Pathogen"] = [pathogen] + [""] * (df.shape[0] - 1)
        df["Genetic Group Type"] = [genetic_group_type] + [""] * (df.shape[0] - 1)
        tables.append(df)

    df = pd.concat(tables, ignore_index=True)

    embedding_name_by_abbreviation = {
        "pca": "PCA",
        "mds": "MDS",
        "t-sne": "t-SNE",
        "umap": "UMAP",
    }
    df["method"] = df["method"].map(embedding_name_by_abbreviation)
    df = df.rename(columns={
        "method": "Method",
        "distance_threshold": "Threshold",
        "normalized_vi": "VI",
    })

    if args.output_csv:
        df.to_csv(args.output_csv, sep=args.separator, header=True)

    if args.output_table:
        df.reset_index().to_latex(
            args.output_table,
            bold_rows=True,
            index=False,
            columns=["Pathogen", "Genetic Group Type", "Method", "VI", "Threshold"],
        )
