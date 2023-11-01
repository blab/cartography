"""Concatenate two or more tables as data frames.
"""
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", nargs="+", help="tables to concatenate")
    parser.add_argument("--separator", default="\t", help="separator between columns in the given tables")
    parser.add_argument("--sort-by", help="what to sort the dataframe by (optional)")
    parser.add_argument("--dataset-name", help="name of dataset to add as column")
    parser.add_argument("--output", help="concatenated table")
    parser.add_argument("--output-latex", help="concatenated table in LaTeX format")
    parser.add_argument("--latex-precision", type=int, default=4, help="precision to use for LaTeX table output")

    args = parser.parse_args()

    # Concatenate tables.
    df = pd.concat([
        pd.read_csv(table_file, sep=args.separator)
        for table_file in args.tables
    ], ignore_index=True)

    if args.sort_by is not None:
        df = df.sort_values(by=[args.sort_by])
        cols_to_order = [args.sort_by]
        new_columns = cols_to_order + (df.columns.drop(cols_to_order).tolist())
        df = df[new_columns]

    if args.dataset_name is not None:
        df["dataset_name"] = args.dataset_name

    df.to_csv(args.output, sep=args.separator, header=True, index=False)

    if args.output_latex:
        df.to_latex(
            args.output_latex,
            bold_rows=True,
            index=False,
            float_format=f"%.{args.latex_precision}f",
        )
