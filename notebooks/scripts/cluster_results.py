import argparse
from augur.utils import write_json, read_node_data
import Bio.SeqIO
from collections import OrderedDict
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

from Helpers import get_euclidean_data_frame, get_PCA_feature_matrix

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
        #threshold_df = pd.read_csv(args.threshold_information)    threshold_df.loc[threshold_df['embedding'] == args.method][args.column_threshold].values.tolist()[0]
        distance_thresholds = args.threshold_information 
    else:
        distance_thresholds = np.arange(0,20,2)

    default_tuned_values = []
    list_of_embedding_strings = ["t-sne", "umap", "mds", "pca"] #["t-SNE","UMAP","MDS", "PCA"]
    embedding_class = [TSNE, UMAP, MDS, PCA]

    # reading in the distance matrix and node data

    distance_matrix = pd.read_csv(args.distance_matrix, index_col=0)
    sequence_names = distance_matrix.index.values.tolist()


    # parameters for the methods taken from the exhaustive grid search 
    embedding_parameters = {
        "metric": "precomputed",
        "perplexity": 30,
        "learning_rate": 500.0,
        "square_distances" : True
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "init": "spectral",
        'n_neighbors': 25, 
        "min_dist": 0.05
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "dissimilarity": "precomputed",
        "n_components" : 2,
        "n_init" : 2
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "n_components" : 10,
        "svd_solver" : "full"
    }
    default_tuned_values.append(embedding_parameters)
    
    
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

    numbers = get_PCA_feature_matrix(args.alignment, sequence_names)

# k fold analysis with thresholds
    random_state = 12883823
    rkf = RepeatedKFold(n_splits=2,  n_repeats=3, random_state=random_state)
    k = 0
    total_list_methods = []
    for training_index, validation_index in rkf.split(clade_annotations["strain"].values.tolist()): 
        print(len(training_index))
        i = 0
        print("here " + str(k))
        for embed in default_tuned_values:
            print(i)
            for distance_threshold in distance_thresholds:

                if(list_of_embedding_strings[i] == "pca"):
                    #performing PCA on my pandas dataframe
                    numbers_subset = numbers[training_index]
                    pca = PCA(**embed) #can specify n, since with no prior knowledge, I use None
                    training_embedding = pca.fit_transform(numbers_subset)
                else:
                    # Subset distance matrix to training indices.
                    training_distance_matrix = distance_matrix[training_index][:, training_index]

                    # Embed training distance matrix.
                    embedder = embedding_class[i](**embed)
                    training_embedding = embedder.fit_transform(training_distance_matrix)

                list_columns_val = _get_embedding_columns_by_method(list_of_embedding_strings[i])
                val_df = pd.DataFrame(training_embedding, columns=list_columns_val)
                val_df[["strain", "clade_membership"]] = pd.DataFrame(clade_annotations[["strain", "clade_membership"]].values[training_index].tolist())

                KDE_df_normal = get_euclidean_data_frame(sampled_df=val_df, column_for_analysis="clade_membership", embedding="method", column_list=_get_embedding_columns_by_method(list_of_embedding_strings[i]))
                
                
                distance_threshold = float(distance_threshold)
                
                clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=distance_threshold)
                clusterer.fit(val_df[_get_embedding_columns_by_method(list_of_embedding_strings[i])])
                val_df[f"{list_of_embedding_strings[i]}_label_{k}"] = clusterer.labels_.astype(str)
                list_columns = _get_embedding_columns_by_method(list_of_embedding_strings[i])
                list_columns.extend(["strain", f"{list_of_embedding_strings[i]}_label_{k}"])
                
                KDE_df_cluster = get_euclidean_data_frame(sampled_df=val_df[list_columns], column_for_analysis=f"{list_of_embedding_strings[i]}_label_{k}", embedding=list_of_embedding_strings[i], column_list=_get_embedding_columns_by_method(list_of_embedding_strings[i]))

                confusion_matrix_val = confusion_matrix(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
                matthews_cc_val = matthews_corrcoef(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
                method_dict = {}

                method_dict["method"] = list_of_embedding_strings[i]
                method_dict["distance_threshold_number"] = f"{list_of_embedding_strings[i]}_label_{k}"
                method_dict["confusion_matrix_training"] = confusion_matrix_val
                method_dict["matthews_cc_training"] = matthews_cc_val
                method_dict["num_undefined_training"] = (val_df[f"{list_of_embedding_strings[i]}_label_{k}"].values == '-1').sum()
                

                if(list_of_embedding_strings[i] == "pca"):
                    #performing PCA on my pandas dataframe
                    numbers_subset = numbers[validation_index]
                    pca = PCA(**embed) #can specify n, since with no prior knowledge, I use None
                    validation_embedding = pca.fit_transform(numbers_subset)
                else:
                    # Subset distance matrix to validation indices.
                    validation_distance_matrix = distance_matrix[validation_index][:, validation_index]

                    # Embed validation distance matrix.
                    validation_embedding = embedder.fit_transform(validation_distance_matrix)

                val_df = pd.DataFrame(validation_embedding, columns=list_columns_val)
                val_df[["strain", "clade_membership"]] = pd.DataFrame(clade_annotations[["strain", "clade_membership"]].values[validation_index].tolist())

                KDE_df_normal = get_euclidean_data_frame(sampled_df=val_df, column_for_analysis="clade_membership", embedding="method", column_list=_get_embedding_columns_by_method(list_of_embedding_strings[i]))
                
                distance_threshold = float(distance_threshold)
                
                clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=distance_threshold)
                clusterer.fit(val_df[_get_embedding_columns_by_method(list_of_embedding_strings[i])])
                val_df[f"{list_of_embedding_strings[i]}_label_{k}"] = clusterer.labels_.astype(str)
                list_columns = _get_embedding_columns_by_method(list_of_embedding_strings[i])
                list_columns.extend(["strain", f"{list_of_embedding_strings[i]}_label_{k}"])
                
                KDE_df_cluster = get_euclidean_data_frame(sampled_df=val_df[list_columns], column_for_analysis=f"{list_of_embedding_strings[i]}_label_{k}", embedding=list_of_embedding_strings[i], column_list=_get_embedding_columns_by_method(list_of_embedding_strings[i]))

                confusion_matrix_val = confusion_matrix(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
                matthews_cc_val = matthews_corrcoef(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])


                method_dict["confusion_matrix_validation"] = confusion_matrix_val
                method_dict["matthews_cc_validation"] = matthews_cc_val
                method_dict["threshold"] = distance_threshold
                method_dict["num_undefined_validation"] = (val_df[f"{list_of_embedding_strings[i]}_label_{k}"].values == '-1').sum()



                print(method_dict)
                total_list_methods.append(method_dict)
                k = k + 1
                
            i = i + 1

        full_output_df = pd.DataFrame(total_list_methods)
            
        if(args.output_full):
            full_output_df.to_csv(args.output_full)

        if args.output_figure:
            sns.relplot(data=full_output_df, x="threshold", y="matthews_cc_validation", col="method", kind="scatter")
            #seaborn, tight layout (more whitespace), larger font, plot mean value as sep mark (compare across thresholds and across method)
            plt.savefig(args.output_figure)

        if args.output:
            max_values = []
            for method in list_of_embedding_strings:
                method_dict = dict(full_output_df.groupby("method").get_group(method).iloc[full_output_df.groupby("method").get_group(method).groupby("threshold")["matthews_cc_validation"].mean().argmax()])
                max_values.append(method_dict)
                
            max_df = pd.DataFrame(max_values)
            max_df.to_csv(args.output)
            print(max_df)
            max_thresholds = max_df["threshold"].values.tolist()
        




