import argparse
import pandas as pd
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", nargs="+", help="KDE plots")
    parser.add_argument("--output", help="comma or tab-delimited file of attributes per node of the given tree")
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
    
    df.to_csv(args.output, sep=args.separator, header=True, index=True, columns=["embedding", "MCC", "accuracy_confusion_matrix", "TN", "FN", "TP", "FP", "classifier_threshold", "median_within", "median_between"])