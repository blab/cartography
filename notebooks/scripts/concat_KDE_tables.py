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

    args = parser.parse_args()

    # Concatenate tables.
    df = pd.concat([
        pd.read_csv(table_file, sep=args.separator)
        for table_file in args.tables
    ], ignore_index=True)

    table_index = []
    for i in range(0, len(args.tables)):
        table_index.append(args.disease_names[i])
        table_index.append([" "] * (len(pd.read_csv(args.tables[i], sep=args.separator).index) -1) )

    def flatten(A):
        rt = []
        for i in A:
            if isinstance(i,list): rt.extend(flatten(i))
            else: rt.append(i)
        return rt  

    df.index = flatten(table_index)

    if args.output_csv:
        df.to_csv(args.output_csv, sep=args.separator, header=True, index=True, columns=["embedding", "MCC", "accuracy", "TN", "FN", "TP", "FP", "threshold", "median_within", "median_between"])
    if args.output_table:
        with open(args.output_table, 'w') as f:
            f.write(df.to_markdown(index=True, tablefmt="grid"))