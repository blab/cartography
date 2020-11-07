import sys
sys.path.append("../")

import argparse
from augur.utils import read_node_data
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import sys

from Helpers import get_euclidean_data_frame

if __name__ == "__main__":

    # Initialize parsers
    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--embedding", help="the path to a dataframe csv file OR distance matrix csv for a genetic KDE plot")
    parser.add_argument("--clades", help="a path to the clade status of the different strains in the build")
    parser.add_argument("--metadata", help="a path to a tsv file that contains information about the differentiator column")
    parser.add_argument("--method", required=True, choices = ["pca", "mds", "t-sne", "umap", "genetic"], help="the embedding used")
    parser.add_argument("--embedding-columns", nargs="+", required=True, help="list of the columns to use as coordinates from the embedding data frame")
    parser.add_argument("--differentiator-column", default="clade_membership", help="string name of the column to differentiate by (clade, host, etc)")
    parser.add_argument("--output-figure", help="path for outputting as a PNG")
    parser.add_argument("--output-dataframe", help="path for outputting as a dataframe")
    parser.add_argument("--output-metadata", help="return the Matthews Correlation Coefficient, Median within and between thresholds and accuracy values for the KDE density plot")

    args = parser.parse_args()

    #Error Handling
    if args.output_figure is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)

    if args.clades is None and args.method is None:
        print("You must specify one of the two", file=sys.stderr)
        sys.exit(1)

    #read in embedding file
    embedding_df = pd.read_csv(args.embedding, index_col=0)

    #read in node data
    if args.clades is not None:
        node_data = read_node_data(args.clades)

        sequences_by_name = list(embedding_df.index)

        sequences_list = []
        for sequence in sequences_by_name:
            if sequence in node_data["nodes"]:
                sequences_list.append(sequence)

        # Build a data frame of clade annotations per strain in the same order
        # as the sequences and the subsequent distance matrix.
        clade_annotations = pd.DataFrame([
            {"strain": sequence_name, args.differentiator_column: node_data["nodes"][sequence_name][args.differentiator_column]}
            for sequence_name in sequences_list
        ])
    elif args.metadata is not None:
        node_data = pd.read_csv(args.metadata, sep='\t')
        node_data.index = node_data["strain"]
        node_data.drop("strain", axis=1)
        node_dict = node_data.transpose().to_dict()

        clade_annotations = pd.DataFrame([
            {"strain": sequence, args.differentiator_column: node_dict[sequence][args.differentiator_column]}
            for sequence in node_data.index
        ])

    #initializing scaler
    scaler = StandardScaler()

    if args.method != "genetic":
        merged_df = embedding_df.merge(clade_annotations, on="strain")
        KDE_df = get_euclidean_data_frame(sampled_df=merged_df, column_list=args.embedding_columns, column_for_analysis=args.differentiator_column, embedding=args.method)
        KDE_df["scaled_distance"] = scaler.fit_transform(KDE_df["distance"].values.reshape(-1, 1)).flatten()

    if args.method == "genetic":
        embedding_df.columns = embedding_df.index
        indices_to_drop = embedding_df[~embedding_df.index.isin(clade_annotations["strain"])].dropna(how = 'all')
        embedding_df = embedding_df[embedding_df.index.isin(clade_annotations["strain"])].dropna(how = 'all')
        embedding_df = embedding_df.drop(indices_to_drop.index, axis=1)
        embedding_df["strain"] = embedding_df.index

        merged_df = embedding_df.merge(clade_annotations, on="strain")

        KDE_df = get_euclidean_data_frame(sampled_df=merged_df, column_for_analysis=args.differentiator_column, embedding=args.method)
        KDE_df["scaled_distance"] = scaler.fit_transform(np.array(KDE_df["distance"]).reshape(-1, 1)).flatten()

    # Use a support vector machine classifier to identify an optimal threshold
    # to distinguish between within and between class pairs. Assign higher
    # weight to "within group" data, to address unbalanced nature of these data
    # with respect to the "between group" data.
    #
    # See also: https://scikit-learn.org/stable/modules/svm.html#svm-classification
    classifier = LinearSVC(
        dual=False,
        random_state=0,
        class_weight={1: 5},
        verbose=0
    )
    X = np.array(KDE_df["scaled_distance"]).reshape(-1,1)
    y = KDE_df["clade_status"].astype(float).values
    classifier.fit(X, y)

    # Find the SVM's threshold between the two given classes by passing the
    # range of possible scaled distance values (effectively z-scores) to the
    # classifier's decision function. When the resulting "scores" are >0, the
    # positive class (e.g., "within clade") are predicted for the corresponding
    # input values. Since negative scaled distances indicate samples that are
    # closer together and the class labels is 1 for "within clade" samples, we
    # look for the highest valued scaled distance for which the decision
    # function returns a positive value. See the documentation for more:
    # https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html#sklearn.svm.LinearSVC.decision_function
    x_range = np.linspace(-3, 3, 1000)
    z = classifier.decision_function(x_range.reshape(-1, 1))
    try:
        classifier_threshold = x_range[np.argwhere(z > 0)[-1]][0]
    except:
        classifier_threshold = "NaN"

    # Estimate group labels using the same input data used to train the SVM.
    estimated_clade_status = classifier.predict(X)

    print(KDE_df["clade_status"].value_counts().sort_values(ascending=False))
    print(list(set(estimated_clade_status)))

    #creating metrics for quantifying patterns within the graph
    confusion_matrix_val = confusion_matrix(y, estimated_clade_status)
    print(confusion_matrix_val)

    accuracy = classifier.score(X, y)
    matthews_cc_val = matthews_corrcoef(y, estimated_clade_status)
    median_within = np.median(KDE_df.query("clade_status == 1")["scaled_distance"])
    median_between = np.median(KDE_df.query("clade_status == 0")["scaled_distance"])

    #create metadata dataframe
    if args.output_metadata is not None:
        metadata_df = pd.DataFrame([[matthews_cc_val, accuracy, median_within, median_between, classifier_threshold, args.method, confusion_matrix_val[0][0], confusion_matrix_val[1][0], confusion_matrix_val[1][1], confusion_matrix_val[0][1]]], columns=["MCC", "accuracy", "median_within", "median_between", "threshold", "embedding", "TN", "FN", "TP", "FP"]).round(3)
        metadata_df.to_csv(args.output_metadata, index=False)

    if args.output_dataframe is not None:
        KDE_df.to_csv(args.output_dataframe)

    if args.output_figure is not None:

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax = sns.kdeplot(KDE_df.query("clade_status == 1")["scaled_distance"], label="Same clade", ax=ax)
        ax = sns.kdeplot(KDE_df.query("clade_status == 0")["scaled_distance"], label="Different clade", ax=ax)

        ax.axvline(x=classifier_threshold, label="SVC threshold", color="#000000", alpha=0.5)

        ax.set_xlabel("Scaled Euclidean distance from embedding")
        ax.set_ylabel("KDE density")
        ax.legend(frameon=False)

        fig.suptitle(args.method + ' KDE Plot', fontsize=16)
        sns.despine()

        plt.tight_layout()
        plt.savefig(args.output_figure)
