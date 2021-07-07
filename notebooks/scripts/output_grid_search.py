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

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--methods", nargs="+", help="methods in the file")
    parser.add_argument("--grid-search-total", help="TSV file with the grid search data")
    parser.add_argument("--output", nargs=2, help="the path where the best thresholds will be saved.")
    parser.add_argument("--output-hyperparameters", nargs=2, help="the path where the best parameters will be saved. ")
    parser.add_argument("--output-figure-HDBSCAN", help="PNG with the results displayed graphically for HDBSCAN thresholds")
    parser.add_argument("--output-figure-grid-search", help="PNG with the results displayed graphically for grid search")

    args = parser.parse_args()

    df = pd.read_csv(args.grid_search_total, sep="\t")

    if args.output_figure_HDBSCAN:
        #TODO: filter dataframe to best set of parameters for t-sne and umap 
        grouped_df = df.groupby(["method", "distance_threshold"])
        maximums = grouped_df.max()
        maximums = maximums.reset_index() 
        sns.relplot(data=maximums, x="distance_threshold", y="validation_mcc", col="method", kind="scatter")
        plt.savefig(args.output_figure_HDBSCAN)
        
    if args.output_figure_grid_search is not None:

        sns.set_theme()

        fig = plt.figure(figsize=(16, 8), constrained_layout=False)
        gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.4, wspace=0.6)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])

        # Creates two subplots and unpacks the output array immediately
        sns.scatterplot(x='learning_rate', y='training_mcc', data=df, hue='perplexity', palette="Set1", ax=ax1)
        ax1.set_xlabel("Learning Rate")
        ax1.set_ylabel("MCC")
        ax1.set_title('TSNE')

        sns.scatterplot(x='perplexity', y='training_mcc', data=df, hue='learning_rate', palette="Set1", ax=ax2)
        ax2.set_xlabel("Perplexity")
        ax2.set_ylabel("MCC")
        ax2.set_title('TSNE')

        sns.scatterplot(x='n_neighbors', y='training_mcc', data=df, hue='min_dist', palette="Set1", ax=ax3)
        ax3.set_xlabel("N Neighbors")
        ax3.set_ylabel("MCC")
        ax3.set_title("UMAP")

        sns.scatterplot(x='min_dist', y='training_mcc', data=df, hue='n_neighbors', palette="Set1", ax=ax4)
        ax4.set_xlabel("Minimum Distance")
        ax4.set_ylabel("MCC")
        ax4.set_title("UMAP")


        ax1.set_ylim(0,1)
        ax2.set_ylim(0,1)
        ax3.set_ylim(0,1)
        ax4.set_ylim(0,1)

        plt.savefig(args.output_figure_grid_search)

    if args.output is not None:
        #make this a dataframe
        max_values = []
        for method in args.methods:
            method_dict = dict(df.groupby("method").get_group(method).iloc[df.groupby("method").get_group(method).groupby("distance_threshold")["validation_mcc"].mean().argmax()])
            max_values.append(method_dict)
                
        max_df = pd.DataFrame(max_values)
        max_index = max_df["method"].values.tolist()
        max_thresholds = max_df["distance_threshold"].values.tolist()
        
        max_df.to_csv(args.output[0])

        df_TSNE = df[df.method == 't-sne'].dropna(axis = 1)
        df_UMAP = df[df.method == 'umap'].dropna(axis = 1)

        TSNE_grouped = pd.DataFrame(df_TSNE.groupby(["perplexity", "learning_rate"])['training_mcc'].mean())
        tsne_val = TSNE_grouped.iloc[TSNE_grouped["training_mcc"].argmax()]

        UMAP_grouped = pd.DataFrame(df_UMAP.groupby(["n_neighbors", "min_dist"])['training_mcc'].mean())
        umap_val = UMAP_grouped.iloc[UMAP_grouped["training_mcc"].argmax()]

        file = open(args.output[1], "w")

        file.write("tsne perplexity: " + str(tsne_val.name[0]) + "\n" + "tsne learning_rate: " + str(tsne_val.name[1]) + "\n" + "mcc best value: " + str(tsne_val.values[0]) + "\n")

        file.write("umap nearest_neighbors: " + str(umap_val.name[0]) + "\n" + "umap min_dist: " + str(umap_val.name[1]) + "\n" + "mcc best value: " + str(umap_val.values[0]))

        file.write("\n".join([str(max_index[i]) + " best threshold is " + str(max_thresholds[i]) + "\n" for i in range(0,len(max_thresholds))]))

        file.close()

    if args.output_hyperparameters is not None:

        max_values = []
        for method in args.methods:
            method_dict = dict(df.groupby("method").get_group(method).iloc[df.groupby("method").get_group(method).groupby("distance_threshold")["validation_mcc"].mean().argmax()])
            max_values.append(method_dict)
                
        max_df = pd.DataFrame(max_values)
        max_index = max_df["method"].values.tolist()
        max_thresholds = max_df["distance_threshold"].values.tolist()
        
        max_df.to_csv(args.output_hyperparameters[0])

        df_TSNE = df[df.method == 't-sne'].dropna(axis = 1)
        df_UMAP = df[df.method == 'umap'].dropna(axis = 1)

        TSNE_grouped = pd.DataFrame(df_TSNE.groupby(["perplexity", "learning_rate"])['training_mcc'].mean())
        tsne_val = TSNE_grouped.iloc[TSNE_grouped["training_mcc"].argmax()]

        UMAP_grouped = pd.DataFrame(df_UMAP.groupby(["n_neighbors", "min_dist"])['training_mcc'].mean())
        umap_val = UMAP_grouped.iloc[UMAP_grouped["training_mcc"].argmax()]

        file = open(args.output_hyperparameters[1], "w")

        file.write("tsne perplexity: " + str(tsne_val.name[0]) + "\n" + "tsne learning_rate: " + str(tsne_val.name[1]) + "\n" + "mcc best value: " + str(tsne_val.values[0]) + "\n")

        file.write("umap nearest_neighbors: " + str(umap_val.name[0]) + "\n" + "umap min_dist: " + str(umap_val.name[1]) + "\n" + "mcc best value: " + str(umap_val.values[0]))

        file.write("\n".join([str(max_index[i]) + " best threshold is " + str(max_thresholds[i]) + "\n" for i in range(0,len(max_thresholds))]))

        file.close()

