import argparse
from augur.utils import write_json, read_node_data
import hdbscan
import matplotlib.pyplot as plt
import math
import matplotlib as mpl
import numpy as np
import pandas as pd
import re
from scipy.spatial.distance import pdist
import seaborn as sns
from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from sklearn.model_selection import RepeatedKFold
import sys
from umap import UMAP

import sys
sys.path.append("../notebooks/scripts/")

from Helpers import get_PCA_feature_matrix

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--distance-matrix", help="csv file with the distance matrix")
    parser.add_argument("--alignment", help="FASTA file with the alignment")
    parser.add_argument("--clades", help="json file containing information about clade membership")
    parser.add_argument("--column-metadata", default="clade_membership", help="the column which contains the clade information")
    parser.add_argument("--threshold-information", nargs="+", help="the distance threshold values to be used on HDBSCAN. if not provided, it will run without.")
    parser.add_argument("--output", help="the csv path where the best label of the data and strain name per method will be saved.")
    parser.add_argument("--output-full", help="the csv path where the full list of accuracy data per threshold will be saved")
    parser.add_argument("--output-figure", help="PNG with the MCC values displayed graphically per method")


    args = parser.parse_args()

    if args.output is None and args.output_cluster_labels is not None:
        print("You must select output in order to get cluster labels", file=sys.stderr)
        sys.exit(1)


    def _get_embedding_columns_by_method(method):
        if method in ("pca"):
            return list(f"{method}1 {method}2 {method}3 {method}4 {method}5 {method}6 {method}7 {method}8 {method}9 {method}10".split())
        if method in ("mds"):
            return list(f"{method}1 {method}2".split())
        if method in ("t-sne"):
            return list("tsne_x tsne_y".split())
        else:
            return list(f"{method}_x {method}_y".split())

    if(args.threshold_information is not None):
        distance_thresholds = args.threshold_information 
    else:
        distance_thresholds = np.arange(0,20,2)

    default_tuned_values = []
    list_of_embedding_strings = "MDS" #["t-SNE","UMAP","MDS", "PCA"]
    embedding_class = MDS

    # reading in the distance matrix and node data

    distance_matrix = pd.read_csv(args.distance_matrix, index_col=0)
    sequence_names = distance_matrix.index.values.tolist()


    # parameters for the methods taken from the exhaustive grid search 

    embedding_parameters = {
        "dissimilarity": "precomputed",
        "n_components" : 2,
        "n_init" : 2,
        "n_jobs" : 1
    }

    # creating dataframe of clade information 
    
    node_data = read_node_data(args.clades)
    clade_annotations = pd.DataFrame([
    {"strain": sequence_name, "clade_membership": node_data["nodes"][sequence_name][args.column_metadata]}
    for sequence_name in sequence_names if sequence_name in node_data["nodes"]
    ])
    
    strains_df = pd.DataFrame(distance_matrix.index.values.tolist(), columns=["strain"])
    clade_annotations = clade_annotations.merge(strains_df, on="strain")

    distance_matrix.columns = distance_matrix.index
    indices_to_drop = distance_matrix[~distance_matrix.index.isin(clade_annotations["strain"])]
    distance_matrix = distance_matrix[distance_matrix.index.isin(clade_annotations["strain"])].dropna(how = 'all')
    distance_matrix = distance_matrix.drop(indices_to_drop.index, axis=1)
    distance_matrix = distance_matrix.to_numpy()
    

# k fold analysis with thresholds

#find main cluster (not_outliers), the rest are outliers
    random_state = 12883823
    rkf = RepeatedKFold(n_splits=2,  n_repeats=3, random_state=random_state)
    k = 0
    total_list_methods = []
    for training_index, validation_index in rkf.split(clade_annotations["strain"].values.tolist()): 
        print("here " + str(k))
        i=0
        for distance_threshold in distance_thresholds:

            # Subset distance matrix to training indices.
            training_distance_matrix = distance_matrix[training_index][:, training_index]

            training_annotations = clade_annotations.iloc[training_index,:]
            # Embed training distance matrix.
            embedder = embedding_class(**embedding_parameters)
            training_embedding = embedder.fit_transform(training_distance_matrix)

            list_columns_val = _get_embedding_columns_by_method(list_of_embedding_strings)
            val_df = pd.DataFrame(training_embedding, columns=list_columns_val)
            val_df[["strain", "clade_membership"]] = pd.DataFrame(clade_annotations[["strain", "clade_membership"]].values[training_index].tolist())

            distance_threshold = float(distance_threshold)
            
            clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=distance_threshold)
            clusterer.fit(val_df[_get_embedding_columns_by_method(list_of_embedding_strings)])
            val_df[f"{list_of_embedding_strings}_label_{k}"] = clusterer.labels_.astype(str)
            list_columns = _get_embedding_columns_by_method(list_of_embedding_strings)
            list_columns.extend(["strain", f"{list_of_embedding_strings}_label_{k}"])

            index_max = val_df[f"{list_of_embedding_strings}_label_{k}"].value_counts().idxmax()
            
            val_df["outlier_status_predicted"] = val_df[f"{list_of_embedding_strings}_label_{k}"].apply(lambda label: "outlier" if label!=index_max else "not_outlier")
            

            confusion_matrix_val = confusion_matrix(training_annotations["clade_membership"], val_df["outlier_status_predicted"])
            matthews_cc_val = matthews_corrcoef(training_annotations["clade_membership"], val_df["outlier_status_predicted"])
            method_dict = {}

            method_dict["method"] = list_of_embedding_strings
            method_dict["confusion_matrix_training"] = confusion_matrix_val
            method_dict["matthews_cc_training"] = matthews_cc_val
            method_dict["num_undefined_training"] = (val_df[f"{list_of_embedding_strings}_label_{k}"].values == '-1').sum()
            
            validation_distance_matrix = distance_matrix[validation_index][:, validation_index]

            validation_annotations = clade_annotations.iloc[validation_index,:]
            # Embed validation distance matrix.
            embedder = embedding_class(**embedding_parameters)
            validation_embedding = embedder.fit_transform(validation_distance_matrix)

            list_columns_val = _get_embedding_columns_by_method(list_of_embedding_strings)
            val_df = pd.DataFrame(validation_embedding, columns=list_columns_val)
            val_df[["strain", "clade_membership"]] = pd.DataFrame(clade_annotations[["strain", "clade_membership"]].values[validation_index].tolist())

            distance_threshold = float(distance_threshold)
            
            clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=distance_threshold)
            clusterer.fit(val_df[_get_embedding_columns_by_method(list_of_embedding_strings)])
            val_df[f"{list_of_embedding_strings}_label_{k}"] = clusterer.labels_.astype(str)
            list_columns = _get_embedding_columns_by_method(list_of_embedding_strings)
            list_columns.extend(["strain", f"{list_of_embedding_strings}_label_{k}"])
            
            index_max = val_df[f"{list_of_embedding_strings}_label_{k}"].value_counts().idxmax()
            
            val_df["outlier_status_predicted"] = val_df[f"{list_of_embedding_strings}_label_{k}"].apply(lambda label: "outlier" if label!=index_max else "not_outlier")
            

            confusion_matrix_val = confusion_matrix(validation_annotations["clade_membership"], val_df["outlier_status_predicted"])
            matthews_cc_val = matthews_corrcoef(validation_annotations["clade_membership"], val_df["outlier_status_predicted"])
            method_dict = {}

            method_dict["method"] = list_of_embedding_strings
            method_dict["distance_threshold_number"] = f"{list_of_embedding_strings}_label_{k}"
            method_dict["distance_threshold"] = distance_threshold
            method_dict["confusion_matrix_validation"] = confusion_matrix_val
            method_dict["matthews_cc_validation"] = matthews_cc_val
            method_dict["num_undefined_validation"] = (val_df[f"{list_of_embedding_strings}_label_{k}"].values == '-1').sum()
            
            print(method_dict)
            total_list_methods.append(method_dict)
            i=i+1
        k = k + 1

    full_output_df = pd.DataFrame(total_list_methods)
            
    if(args.output_full):
        full_output_df.to_csv(args.output_full)

    if args.output_figure:
        sns.relplot(data=full_output_df, x="distance_threshold", y="matthews_cc_validation", col="method", kind="scatter")
        plt.savefig(args.output_figure)

    if args.output:
        max_values = []
        for method in list_of_embedding_strings:
            # method_val: groups the full dataframe by method, returns the dataframes grouped by method
            method_val = full_output_df.groupby("method").get_group(method)
            method_dict = dict(method_val.iloc[method_val.groupby("distance_threshold")["matthews_cc_validation"].mean().argmax()])
            max_values.append(method_dict)
            
        max_df = pd.DataFrame(max_values)
        max_df.to_csv(args.output)
        print(max_df)
        max_thresholds = max_df["distance_threshold"].values.tolist()
        




