#!/usr/bin/env python3
import argparse
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    with open(args.json, "r", encoding="utf-8") as fh:
        schema = json.load(fh)

    columns_to_drop = [
        "parent_name",
        "branch_length",
        "pca_label",
        "mds_label",
        "t-sne_label",
        "umap_label",
        "genetic_label",
        "mds3",
    ]

    for dataset in schema["data"]:
        if "values" in dataset:
            for record in dataset["values"]:
                for column in columns_to_drop:
                    if column in record:
                        del record[column]

                for key, value in record.items():
                    if type(value) is float:
                        record[key] = round(value, 3)

    with open(args.output, "w", encoding="utf-8") as oh:
        json.dump(
            schema,
            oh,
            indent=0,
            separators=(",", ":"),
        )
