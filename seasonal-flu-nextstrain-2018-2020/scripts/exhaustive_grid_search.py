"""Exhastuve grid search for parameters for TSNE and UMAP"""
import argparse
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC
from umap import UMAP

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--distance-matrix", help="csv file with the distance matrix")
    parser.add_argument("--node-data", help="csv file with the clade_membership - that MUST be the name of the column.")
    parser.add_argument("--n-neighbors", nargs="+", type=int, help="list of values that the search should use")
    parser.add_argument("--min-dist", nargs="+", type=float, help="list of values that the search should use")
    parser.add_argument("--perplexity", nargs="+", type=float, help="list of values that the search should use")
    parser.add_argument("--learning-rate", nargs="+", type=float, help="list of values that the search should use")
    parser.add_argument("--n-repeats", type=int, help="the number of times the k fold generator should repeat the k fold")
    parser.add_argument("--output", help="the path where the best parameters will be saved.")
    parser.add_argument("--output-metadata", help="the path where the grid search data will be saved.")
    parser.add_argument("--output-figure", help="PNG with the results displayed graphically")


    args = parser.parse_args()

    def assign_clade_status_to_pairs(clade_annotations, index):
        """Assign clade status to all pairs in the given list of indices and the given data frame of clade annotations.
        
        Outputs a vector in condensed distance matrix format such that all nonredundant pairs of strains are represented.
        
        """
        clade_statuses = []
        for i in range(len(index)):
            for j in range(i + 1, len(index)):
                same_clade = clade_annotations.loc[index[i], "clade_membership"] == clade_annotations.loc[index[j], "clade_membership"]
                clade_statuses.append(int(same_clade))
                
        return np.array(clade_statuses)

    tuned_parameter_values = []
    list_of_embeddings = [TSNE, UMAP]
    default_tuned_values = []
    list_of_embeddings_strings = ["t-SNE", "UMAP"]

    embedding_parameters = {
        "metric": "precomputed",
    }
    tuned_parameters_TSNE = {
        "perplexity": args.perplexity, #[15, 30, 100],
        "learning_rate": args.learning_rate #[100.0, 200.0, 500.0, 1000.0]
    }

    tuned_parameter_values.append(tuned_parameters_TSNE)
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "init": "spectral",
    }
    tuned_parameters_UMAP = {
        "n_neighbors" : args.n_neighbors, #[25, 100, 200],
        "min_dist" : args.min_dist #[.05, .5]
    }
    tuned_parameter_values.append(tuned_parameters_UMAP)
    default_tuned_values.append(embedding_parameters)

    # reading in the distance matrix and node data
 
    distance_matrix = pd.read_csv(args.distance_matrix, index_col=0)
    distance_matrix.reset_index(drop=True)
    distance_matrix = distance_matrix.to_numpy()

    node_df = pd.read_csv(args.node_data, sep="\t")
    clade_annotations = node_df[["strain", "clade_membership"]]

    sequence_names = node_df["strain"].values.tolist()
    random_state = 12883823
    rkf = RepeatedKFold(n_splits=2, n_repeats=args.n_repeats, random_state=random_state)

    grid_search_results = []
    for training_index, validation_index in rkf.split(sequence_names): 
        i = 0
        scaler = StandardScaler()
        for embed in tuned_parameter_values:
            keys, values = zip(*embed.items())
            experiments = [dict(zip(keys, v)) for v in itertools.product(*values)]
            
            for experiment in experiments:
                method_dict = default_tuned_values[i].copy()
                experiment_tuple = [(k, v) for k, v in experiment.items()]
                method_dict.update(experiment_tuple)

                embedder = list_of_embeddings[i](**method_dict)

                training_distance_matrix = distance_matrix[training_index][:, training_index]
                training_embedding = embedder.fit_transform(training_distance_matrix)

                training_embedding_distances = pdist(training_embedding).reshape(-1, 1)
                training_embedding_distances = scaler.fit_transform(training_embedding_distances).flatten().reshape(-1,1)

                # Assign a binary class to each pair of samples based on their clade memberships.
                # Samples from different clades are assigned 0, samples from the same clade as assigned 1.
                # This vector of binary values will be the output to fit a classifier to.
                # These pairs should be in the same order as the embedding distances above.
                training_clade_status_for_pairs = assign_clade_status_to_pairs(
                    clade_annotations,
                    training_index
                )
                # Use a support vector machine classifier to identify an optimal threshold
                # to distinguish between within and between class pairs.
                # See also: https://scikit-learn.org/stable/modules/svm.html#svm-classification
                classifier = LinearSVC(
                    dual=False,
                    random_state=0,
                    class_weight={1: 5},
                    verbose=0
                )
                X = np.array(training_embedding_distances).reshape(-1,1)
                y = training_clade_status_for_pairs
                classifier.fit(X, y)

                x_range = np.linspace(-3, 3, 1000)
                z = classifier.decision_function(x_range.reshape(-1, 1))
                classifier_threshold = x_range[np.argwhere(z > 0)[-1]][0]

                # Use a SVM to identify an optimal threshold for genetic distances.
                genetic_classifier = make_pipeline(
                    StandardScaler(),
                    LinearSVC(
                        dual=False,
                        random_state=0,
                        class_weight={1: 5},
                        verbose=0
                    )
                )
                genetic_classifier.fit(squareform(training_distance_matrix).reshape(-1, 1), training_clade_status_for_pairs)

                # Subset distance matrix to validation indices.
                validation_distance_matrix = distance_matrix[validation_index][:, validation_index]

                # Embed validation distance matrix.
                validation_embedding = embedder.fit_transform(validation_distance_matrix)

                # Calculate Euclidean distance between pairs of samples in the embedding.
                # The output should be a data frame of distances between pairs.
                validation_embedding_distances = pdist(validation_embedding).reshape(-1, 1)
                validation_embedding_distances = scaler.fit_transform(validation_embedding_distances).flatten().reshape(-1,1)

                # Assign a binary class to each pair of samples based on their clade memberships.
                # Samples from different clades are assigned 0, samples from the same clade as assigned 1.
                # This vector of binary values will be the output to fit a classifier to.
                # These pairs should be in the same order as the embedding distances above.
                validation_clade_status_for_pairs = assign_clade_status_to_pairs(
                    clade_annotations,
                    validation_index
                )

                # Predict and score clade status from embedding distances and the trained classifier.
                # The first argument is the set to predict classifier labels for. The second argument
                # is the list of true labels. The return argument is the mean accuracy of the predicted
                # labels.
                # https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html#sklearn.svm.LinearSVC.score

                confusion_matrix_val = confusion_matrix(validation_clade_status_for_pairs, classifier.predict(validation_embedding_distances))
                matthews_cc_val = matthews_corrcoef(validation_clade_status_for_pairs, classifier.predict(validation_embedding_distances))

                accuracy = classifier.score(
                    validation_embedding_distances,
                    validation_clade_status_for_pairs
                )

                genetic_accuracy = genetic_classifier.score(
                    squareform(validation_distance_matrix).reshape(-1, 1),
                    validation_clade_status_for_pairs
                )

                method_dict["method"] = list_of_embeddings_strings[i]

                method_dict["confusion_matrix"] = confusion_matrix_val
                method_dict["matthews_cc"] = matthews_cc_val
                method_dict["threshold"] = classifier_threshold
                method_dict["accuracy"] = accuracy
                print(method_dict)
                grid_search_results.append(method_dict)
            i = i + 1

        
    df = pd.DataFrame(grid_search_results)
    if args.output_metadata is not None:
        df.to_csv(args.output_metadata)

    if args.output_figure is not None:

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
        df_TSNE = df[df.method == 't-SNE'].dropna(axis = 1)
        df_UMAP = df[df.method == 'UMAP'].dropna(axis = 1)

        TSNE_grouped = pd.DataFrame(df_TSNE.groupby(["perplexity", "learning_rate"])['matthews_cc'].mean())
        tsne_val = TSNE_grouped.iloc[TSNE_grouped["matthews_cc"].argmax()]

        UMAP_grouped = pd.DataFrame(df_UMAP.groupby(["n_neighbors", "min_dist"])['matthews_cc'].mean())
        umap_val = UMAP_grouped.iloc[UMAP_grouped["matthews_cc"].argmax()]

        file = open(args.output, "w")

        file.write("tsne perplexity: " + str(tsne_val.name[0]) + "\n" + "tsne learning_rate: " + str(tsne_val.name[1]) + "\n" + "mcc best value: " + str(tsne_val.values[0]) + "\n")

        file.write("umap nearest_neighbors: " + str(umap_val.name[0]) + "\n" + "umap min_dist: " + str(umap_val.name[1]) + "\n" + "mcc best value: " + str(umap_val.values[0]))

        file.close()

