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
    parser.add_argument("--embedding-columns", nargs=2, required=True, help="list of the two columns to use as coordinates from the embedding data frame")
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
        KDE_df = get_euclidean_data_frame(sampled_df=merged_df, column1=args.embedding_columns[0], column2=args.embedding_columns[1], column_for_analysis=args.differentiator_column, embedding=args.method)
        KDE_df["scaled_distance"] = scaler.fit_transform(KDE_df["distance"].values.reshape(-1, 1)).flatten()

    if args.method == "genetic":
        embedding_df.columns = embedding_df.index
        indices_to_drop = embedding_df[~embedding_df.index.isin(clade_annotations["strain"])].dropna(how = 'all')
        embedding_df = embedding_df[embedding_df.index.isin(clade_annotations["strain"])].dropna(how = 'all')
        embedding_df = embedding_df.drop(indices_to_drop.index, axis=1)
        embedding_df["strain"] = embedding_df.index

        merged_df = embedding_df.merge(clade_annotations, on="strain")

        KDE_df = get_euclidean_data_frame(sampled_df=merged_df, column_for_analysis=args.differentiator_column, embedding=args.method)
        KDE_df["scaled_distance"] = scaler.fit_transform(np.array(KDE_df["distance"]).reshape(-1, 1))

    # Use a support vector machine classifier to identify an optimal threshold
    # to distinguish between within and between class pairs.
    # See also: https://scikit-learn.org/stable/modules/svm.html#svm-classification
    classifier = make_pipeline(
        StandardScaler(),
        LinearSVC(random_state=0, tol=1e-5)
    )

    #fitting the classifier with the scaled distance and the within vs between per strain relationship

    classifier.fit(np.array(KDE_df["scaled_distance"]).reshape(-1,1), KDE_df["clade_status"])
    classifier_threshold = (0.5 - classifier.named_steps["linearsvc"].intercept_) / classifier.named_steps["linearsvc"].coef_[0]
    classifier_threshold = classifier_threshold[0]
    print(KDE_df["clade_status"].value_counts().sort_values(ascending=False))
    print(list(set(classifier.predict(np.array(KDE_df["scaled_distance"]).reshape(-1,1)))))
    #creating metrics for quantifying patterns within the graph

    confusion_matrix_val = confusion_matrix(classifier.predict(np.array(KDE_df["scaled_distance"]).reshape(-1,1)), KDE_df["clade_status"])

    confusion_matrix_number = (confusion_matrix_val[0][0] + confusion_matrix_val[1][1]) / float(len(KDE_df))

    matthews_cc_val = matthews_corrcoef(classifier.predict(np.array(KDE_df["scaled_distance"]).reshape(-1,1)), KDE_df["clade_status"])

    median_within = np.median(KDE_df.query("clade_status == 1")["scaled_distance"])
    
    median_between = np.median(KDE_df.query("clade_status == 0")["scaled_distance"])

    #create metadata dataframe
    if args.output_metadata is not None:
        metadata_df = pd.DataFrame([[matthews_cc_val, confusion_matrix_number, median_within, median_between, classifier_threshold]], columns=["matthews_cc", "accuracy_confusion_matrix", "median_within", "median_between", "classifier_threshold"])
        metadata_df.to_csv(args.output_metadata)

    if args.output_dataframe is not None:
        KDE_df.to_csv(args.output_dataframe)

    if args.output_figure is not None:
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax = sns.kdeplot(KDE_df.query("clade_status == 1")["scaled_distance"], label="Same clade", ax=ax)
        ax = sns.kdeplot(KDE_df.query("clade_status == 0")["scaled_distance"], label="Different clade", ax=ax)

        ax.axvline(x=classifier_threshold, label="SVC threshold", color="#000000", alpha=0.5)

        ax.set_xlabel("Scaled Euclidean distance from embedding")
        ax.set_ylabel("KDE density")

        
        fig.suptitle(args.method + ' KDE Plot - ', fontsize=16)
        sns.despine()

        plt.savefig(args.output_figure)