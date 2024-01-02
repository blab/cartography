#!/usr/bin/env python3


"""
        python3 sars-cov-2-nextstrain-2022-2023/scripts/get_distinct_recombinant_lineages.py \
            --lineages {input.lineages} \
            --lineage-counts {input.counts} \
            --min-count {params.min_count} \
            --output-distinct {output.distinct_lineages} \
            --output-table {output.recombinant_and_parental_lineages}

"""
import argparse
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--lineages", required=True, help="official JSON defining Pango lineages")
    parser.add_argument("--lineage-counts", required=True, help="counts per Pango lineage")
    parser.add_argument("--min-count", default=10, type=int, help="minimum records per lineage required for a recombinant lineage and its parental lineages to be included in the output")
    parser.add_argument("--output-distinct", required=True, help="text file with list of distinct recombinant and parental lineages with a header line")
    parser.add_argument("--output-table", required=True, help="CSV table of recombinant and parental lineages")

    args = parser.parse_args()

    with open(args.lineages, "r", encoding="utf-8") as fh:
        lineages = json.load(fh)

    recombinant_lineages = {
        lineage: sorted(set(parents))
        for lineage, parents in lineages.items()
        if isinstance(parents, list)
    }

    # Load lineage counts.
    lineage_counts = pd.read_csv(args.lineage_counts, sep="\t")
    count_by_lineage = lineage_counts.set_index("Nextclade_pango")["count"].to_dict()

    # Filter recombinant lineages to those for which the recombinant
    # and its parents have the minimum number of records available.
    min_count = args.min_count
    covered_recombinant_lineages = {
        lineage: parents
        for lineage, parents in recombinant_lineages.items()
        if (count_by_lineage.get(lineage, 0) >= min_count and
            count_by_lineage.get(parents[0].replace("*", ""), 0) >= min_count and
            count_by_lineage.get(parents[1].replace("*", ""), 0) >= min_count)
    }

    # Create a list of distinct lineages and a table of recombinant and parental
    # lineages.
    distinct_lineages = set()
    lineage_table = []
    for lineage, parents in covered_recombinant_lineages.items():
        distinct_lineages.add(lineage)
        record = {
            "recombinant_X": lineage,
        }

        for label, parent in zip(("A", "B"), parents):
            simple_parent = parent.replace("*", "")
            distinct_lineages.add(simple_parent)
            record[f"parental_{label}"] = simple_parent

        lineage_table.append(record)

    # Save distinct lineages to file.
    distinct_lineages = sorted(distinct_lineages)
    with open(args.output_distinct, "w", encoding="utf-8") as oh:
        print("Nextclade_pango", file=oh)

        for lineage in distinct_lineages:
            print(lineage, file=oh)

    # Save table of recombinant lineages and their parents.
    pd.DataFrame.from_records(lineage_table).to_csv(
        args.output_table,
        header=True,
        index=False,
    )
