import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", help="the metadata TSV file")
    parser.add_argument("--output", help="the output file")

    args = parser.parse_args()

    df = pd.read_table(args.metadata, na_values=["?"])
    #df["low_quality_clock_deviation"] = df["clock_deviation"].abs() > (1.5 * df["clock_deviation"].std())
    #df["low_quality"] = (df.loc[:, ["QC_missing_data", "QC_mixed_sites", "QC_rare_mutations", "QC_snp_clusters"]] != "good").values.sum(axis=1) > 0
    df["low_quality"] = (df.loc[:, ["QC_missing_data", "QC_mixed_sites", "QC_rare_mutations", "QC_snp_clusters"]] != "good").values.sum(axis=1)
    df.to_csv(args.output, sep="\t")