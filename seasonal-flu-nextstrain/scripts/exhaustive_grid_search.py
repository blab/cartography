"""Exhastuve grid search for parameters for TSNE and UMAP"""
import argparse
import itertools
import hdbscan
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC
from umap import UMAP

import sys
sys.path.append("../notebooks/scripts/")

from Helpers import get_PCA_feature_matrix, get_euclidean_data_frame

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--distance-matrix", help="csv file with the distance matrix")
    parser.add_argument("--alignment", help="FASTA file with the alignment")
    parser.add_argument("--node-data", help="csv file with the clade_membership - that MUST be the name of the column.")
    parser.add_argument("--n-neighbors", nargs="+", type=int, help="list of values that the search should use")
    parser.add_argument("--min-dist", nargs="+", type=float, help="list of values that the search should use")
    parser.add_argument("--perplexity", nargs="+", type=float, help="list of values that the search should use")
    parser.add_argument("--threshold-information", nargs="+", help="the distance threshold values to be used on HDBSCAN. if not provided, it will run without.")
    parser.add_argument("--learning-rate", nargs="+", type=float, help="list of values that the search should use")
    parser.add_argument("--n-repeats", type=int, help="the number of times the k fold generator should repeat the k fold")
    parser.add_argument("--output", help="the path where the best thresholds will be saved.")
    parser.add_argument("--output-hyperparameters", help="the path where the best parameters will be saved. ")
    parser.add_argument("--output-metadata", help="the path where the grid search data will be saved.")
    parser.add_argument("--output-figure-HDBSCAN", help="PNG with the results displayed graphically for HDBSCAN thresholds")
    parser.add_argument("--output-figure-grid-search", help="PNG with the results displayed graphically for grid search")

    args = parser.parse_args()

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
    tuned_parameter_values = []

    # reading in the distance matrix and node data

    distance_matrix = pd.read_csv(args.distance_matrix, index_col=0)
    sequence_names = distance_matrix.index.values.tolist()

    
    # parameters for the methods taken from the exhaustive grid search 
    embedding_parameters = {
        "metric": "precomputed",
        "square_distances" : True
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "init": "spectral",
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "dissimilarity": "precomputed",
        "n_components" : 2,
        "n_init" : 2,
        "n_jobs" : 1
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "n_components" : 10,
        "svd_solver" : "full"
    }
    default_tuned_values.append(embedding_parameters)
    
    
    tuned_parameters_TSNE = {
        "perplexity": args.perplexity, #[15, 30, 100],
        "learning_rate": args.learning_rate, #[100.0, 200.0, 500.0, 1000.0],
        "square_distances" : [True]
    }

    tuned_parameter_values.append(tuned_parameters_TSNE)

    tuned_parameters_UMAP = {
        "n_neighbors" : args.n_neighbors, #[25, 100, 200],
        "min_dist" : args.min_dist #[.05, .5]
    }
    tuned_parameter_values.append(tuned_parameters_UMAP)

        
    tuned_parameters_MDS = {
    }

    tuned_parameter_values.append(tuned_parameters_MDS)

    tuned_parameters_PCA = {
    }
    tuned_parameter_values.append(tuned_parameters_PCA)

    # reading in the distance matrix and node data

    distance_matrix = pd.read_csv(args.distance_matrix, index_col=0)
    sequence_names = distance_matrix.index.values.tolist()

    clade_annotations = pd.read_csv(args.node_data, sep="\t")
    clade_annotations = clade_annotations[["strain", "clade_membership"]]
    #clade_annotations = clade_annotations.merge(pd.DataFrame(sequence_names, columns=["strain"]), on="strain")
    
    distance_matrix.columns = distance_matrix.index
    indices_to_drop = distance_matrix[~distance_matrix.index.isin(clade_annotations["strain"])].dropna(how = 'all')
    distance_matrix = distance_matrix[distance_matrix.index.isin(clade_annotations["strain"])].dropna(how = 'all')
    distance_matrix = distance_matrix.drop(indices_to_drop.index, axis=1)
    sequence_names = distance_matrix.index.values.tolist()
    distance_matrix = distance_matrix.to_numpy()
    #sequence_names = clade_annotations["strain"].values.tolist()
    
    numbers = get_PCA_feature_matrix(args.alignment, sequence_names)

    random_state = 12883823
    rkf = RepeatedKFold(n_splits=2, n_repeats=args.n_repeats, random_state=random_state)
    k=0
    total_list_methods = []
    for training_index, validation_index in rkf.split(clade_annotations["strain"].values.tolist()): 
        print(len(training_index))
        print("here " + str(k))
        for embed in tuned_parameter_values:
            keys, values = zip(*embed.items())
            experiments = [dict(zip(keys, v)) for v in itertools.product(*values)]
            i = 0
            for experiment in experiments:
                method_dict = default_tuned_values[i].copy()
                experiment_tuple = [(k, v) for k, v in experiment.items()]
                method_dict.update(experiment_tuple)

                for distance_threshold in distance_thresholds:

                    if(list_of_embedding_strings[i] == "pca"):
                        #performing PCA on my pandas dataframe
                        numbers_subset = numbers[training_index]
                        pca = PCA(**method_dict) #can specify n, since with no prior knowledge, I use None
                        training_embedding = pca.fit_transform(numbers_subset)
                    else:
                        # Subset distance matrix to training indices.
                        training_distance_matrix = distance_matrix[training_index][:, training_index]

                        # Embed training distance matrix.
                        print(method_dict)
                        embedder = embedding_class[i](**method_dict)
                        training_embedding = embedder.fit_transform(training_distance_matrix)

                    list_columns_val = _get_embedding_columns_by_method(list_of_embedding_strings[i])
                    val_df = pd.DataFrame(training_embedding, columns=list_columns_val)
                    val_df[["strain", "clade_membership"]] = pd.DataFrame(clade_annotations[["strain", "clade_membership"]].values[training_index].tolist())

                    KDE_df_normal = get_euclidean_data_frame(sampled_df=val_df, column_for_analysis="clade_membership", embedding="method", column_list=_get_embedding_columns_by_method(list_of_embedding_strings[i]))
                    
                    
                    distance_threshold = float(distance_threshold)
                    
                    clusterer = hdbscan.HDBSCAN(min_cluster_size=15, cluster_selection_epsilon=distance_threshold)
                    clusterer.fit(val_df[_get_embedding_columns_by_method(list_of_embedding_strings[i])])
                    val_df[f"{list_of_embedding_strings[i]}_label_{k}"] = clusterer.labels_.astype(str)
                    list_columns = _get_embedding_columns_by_method(list_of_embedding_strings[i])
                    list_columns.extend(["strain", f"{list_of_embedding_strings[i]}_label_{k}"])
                    
                    KDE_df_cluster = get_euclidean_data_frame(sampled_df=val_df[list_columns], column_for_analysis=f"{list_of_embedding_strings[i]}_label_{k}", embedding=list_of_embedding_strings[i], column_list=_get_embedding_columns_by_method(list_of_embedding_strings[i]))

                    confusion_matrix_val = confusion_matrix(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
                    matthews_cc_val = matthews_corrcoef(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
                    

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
                    
                    clusterer = hdbscan.HDBSCAN(min_cluster_size=15, cluster_selection_epsilon=distance_threshold)
                    clusterer.fit(val_df[_get_embedding_columns_by_method(list_of_embedding_strings[i])])
                    val_df[f"{list_of_embedding_strings[i]}_label_{k}"] = clusterer.labels_.astype(str)
                    list_columns = _get_embedding_columns_by_method(list_of_embedding_strings[i])
                    list_columns.extend(["strain", f"{list_of_embedding_strings[i]}_label_{k}"])
                    
                    KDE_df_cluster = get_euclidean_data_frame(sampled_df=val_df[list_columns], column_for_analysis=f"{list_of_embedding_strings[i]}_label_{k}", embedding=list_of_embedding_strings[i], column_list=_get_embedding_columns_by_method(list_of_embedding_strings[i]))

                    confusion_matrix_val = confusion_matrix(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])
                    matthews_cc_val = matthews_corrcoef(KDE_df_normal["clade_status"], KDE_df_cluster["clade_status"])

                   # method_dict = {}
                    CV_dict = default_tuned_values[i].copy()
                    CV_dict["method"] = list_of_embedding_strings[i]
                    CV_dict["distance_threshold_number"] = f"{list_of_embedding_strings[i]}_label_{k}"
                    CV_dict["distance_threshold"] = distance_threshold
                    CV_dict["confusion_matrix_training"] = confusion_matrix_val
                    CV_dict["matthews_cc_training"] = matthews_cc_val
                    CV_dict["num_undefined_training"] = (val_df[f"{list_of_embedding_strings[i]}_label_{k}"].values == '-1').sum()
                    CV_dict["confusion_matrix"] = confusion_matrix_val
                    CV_dict["matthews_cc"] = matthews_cc_val

                    print(CV_dict)
                    total_list_methods.append(CV_dict)
                i = i + 1
        k = k+1

        
    df = pd.DataFrame(total_list_methods)
    
    if args.output_metadata is not None:
        df.to_csv(args.output_metadata)

    if args.output_figure_HDBSCAN:
        #TODO: filter dataframe to best set of parameters for t-sne and umap 
        sns.relplot(data=df, x="threshold", y="matthews_cc_validation", col="method", kind="scatter")
        plt.savefig(args.output_figure)
        
    if args.output_figure_grid_search is not None:

        sns.set_theme()

        fig = plt.figure(figsize=(16, 8), constrained_layout=False)
        gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.4, wspace=0.6)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])

        # Creates two subplots and unpacks the output array immediately
        sns.scatterplot(x='learning_rate', y='matthews_cc', data=df, hue='perplexity', palette="Set1", ax=ax1)
        ax1.set_xlabel("Learning Rate")
        ax1.set_ylabel("MCC")
        ax1.set_title('TSNE')

        sns.scatterplot(x='perplexity', y='matthews_cc', data=df, hue='learning_rate', palette="Set1", ax=ax2)
        ax2.set_xlabel("Perplexity")
        ax2.set_ylabel("MCC")
        ax2.set_title('TSNE')

        sns.scatterplot(x='n_neighbors', y='matthews_cc', data=df, hue='min_dist', palette="Set1", ax=ax3)
        ax3.set_xlabel("N Neighbors")
        ax3.set_ylabel("MCC")
        ax3.set_title("UMAP")

        sns.scatterplot(x='min_dist', y='matthews_cc', data=df, hue='n_neighbors', palette="Set1", ax=ax4)
        ax4.set_xlabel("Minimum Distance")
        ax4.set_ylabel("MCC")
        ax4.set_title("UMAP")


        ax1.set_ylim(0,1)
        ax2.set_ylim(0,1)
        ax3.set_ylim(0,1)
        ax4.set_ylim(0,1)

        plt.savefig(args.output_figure)

    if args.output is not None:

        max_values = []
        for method in list_of_embedding_strings:
            method_dict = dict(df.groupby("method").get_group(method).iloc[df.groupby("method").get_group(method).groupby("threshold")["matthews_cc_validation"].mean().argmax()])
            max_values.append(method_dict)
                
        max_df = pd.DataFrame(max_values)
        max_index = max_df["method"].values.tolist()
        max_thresholds = max_df["threshold"].values.tolist()
        
        max_df.to_csv(args.output[1])

        df_TSNE = df[df.method == 't-SNE'].dropna(axis = 1)
        df_UMAP = df[df.method == 'UMAP'].dropna(axis = 1)

        TSNE_grouped = pd.DataFrame(df_TSNE.groupby(["perplexity", "learning_rate"])['matthews_cc'].mean())
        tsne_val = TSNE_grouped.iloc[TSNE_grouped["matthews_cc"].argmax()]

        UMAP_grouped = pd.DataFrame(df_UMAP.groupby(["n_neighbors", "min_dist"])['matthews_cc'].mean())
        umap_val = UMAP_grouped.iloc[UMAP_grouped["matthews_cc"].argmax()]

        file = open(args.output[0], "w")

        file.write("tsne perplexity: " + str(tsne_val.name[0]) + "\n" + "tsne learning_rate: " + str(tsne_val.name[1]) + "\n" + "mcc best value: " + str(tsne_val.values[0]) + "\n")

        file.write("umap nearest_neighbors: " + str(umap_val.name[0]) + "\n" + "umap min_dist: " + str(umap_val.name[1]) + "\n" + "mcc best value: " + str(umap_val.values[0]))

        file.write([str(max_index[i]) + " best threshold is " + str(max_thresholds[i]) + "\n" for i in range(0,len(max_thresholds))])

        file.close()

