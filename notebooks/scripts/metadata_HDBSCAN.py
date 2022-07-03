import argparse
import numpy as np
import pandas as pd

from Helpers import variation_of_information


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--metadata", required=True, help="clades information")
    parser.add_argument("--embedding", required=True, help="a csv file with embedding information for method")
    parser.add_argument("--method", required=True, choices = ["pca", "mds", "t-sne", "umap"], help="the embedding used")
    parser.add_argument("--clade-column", default="clade_membership", help="column storing known clade or group membership to use for accuracy calculations")
    parser.add_argument("--missing-data-value", help="string used to represent missing data values that should be dropped from accuracy calculations")
    parser.add_argument("--analysis-name", help="name of analysis to annotate the accuracy values with. Used when outputs will be concatenated downstream across multiple analyses.")
    parser.add_argument("--cluster-data", help="cluster data from embedding and assign labels given via HDBSCAN")
    parser.add_argument("--cluster-threshold", type=float, help="cluster data from embedding and assign labels given via HDBSCAN. Pass in a threshold.")
    parser.add_argument("--output", required=True, help="outputting a csv file of metadata info for HDBSCAN results")

    args = parser.parse_args()

    embedding_df = pd.read_csv(args.embedding)
    metadata_df = pd.read_csv(args.metadata, sep="\t")
    embedding_df = embedding_df.merge(metadata_df[["strain", args.clade_column]], on="strain")

    if args.cluster_data:
        if args.method == "pca":
            n_components = 10
            columns = [
                f"pca{i}"
                for i in range(1, n_components + 1)
            ]
        elif args.method == "mds":
            parameters = pd.read_csv(args.cluster_data)
            n_components = int(parameters["components"].values[0])
            columns = [
                f"mds{i}"
                for i in range(1, n_components + 1)
            ]
        else:
            method = args.method.replace("-", "")
            columns = [f"{method}_x", f"{method}_y"]

        print(f"Calculate VI score with the following columns: {columns}")

    if args.missing_data_value:
        embedding_df[args.clade_column] = embedding_df[args.clade_column].replace(
            args.missing_data_value,
            np.NaN,
        )
        total_rows = embedding_df.shape[0]
        embedding_df.dropna(subset=[args.clade_column], inplace=True)
        non_missing_rows = embedding_df.shape[0]
        print(f"Dropped {total_rows - non_missing_rows} missing values in the clade column '{args.clade_column}'.")

    # Determine clade status from automated cluster labels.
    # KDE_df_cluster = get_euclidean_data_frame(
    #     sampled_df=embedding_df,
    #     column_for_analysis=f"{args.method}_label",
    #     embedding=args.method,
    #     column_list=columns
    # )

    # # Determine clade status from pre-assigned clade or group membership.
    # KDE_df_normal = get_euclidean_data_frame(
    #     sampled_df=embedding_df,
    #     column_for_analysis=args.clade_column,
    #     embedding=args.method,
    #     column_list=columns
    # )

    # Calculate accuracy of automated cluster labels compared to pre-assigned
    # clade or group membership.
    # confusion_matrix_val = confusion_matrix(
    #     KDE_df_normal["clade_status"],
    #     KDE_df_cluster["clade_status"]
    # )
    # matthews_cc_val = matthews_corrcoef(
    #     KDE_df_normal["clade_status"],
    #     KDE_df_cluster["clade_status"]
    # )

    # variation of information score
    cdict = embedding_df[["strain", f"{args.method}_label"]].set_index("strain")
    clade = cdict.groupby([f"{args.method}_label"])
    list_clades = [clade.get_group(x) for x in clade.groups]
    predicted_values = [list(lists.index) for lists in list_clades]

    ldict = embedding_df[["strain", args.clade_column]].set_index("strain")
    clade = ldict.groupby([args.clade_column])
    list_clades = [clade.get_group(x) for x in clade.groups]
    actual_values = [list(lists.index) for lists in list_clades]



    VI_normalized = variation_of_information(predicted_values, actual_values, normalized=True)

    if args.cluster_threshold is not None:
        cluster_threshold = args.cluster_threshold
    elif args.cluster_data is not None:
        max_df = pd.read_csv(args.cluster_data)
        cluster_threshold = float(max_df['distance_threshold'].values.tolist()[0])
    else:
        pass

    output_df = pd.DataFrame(
        [
            [
                args.method,
                VI_normalized,
                cluster_threshold,
            ]
        ],
        columns=["embedding", "normalized_VI", "threshold"]
    ).round(3)

    if args.analysis_name:
        output_df["analysis_name"] = args.analysis_name

    output_df.to_csv(args.output, index=False)
