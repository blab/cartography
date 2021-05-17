"""Exhastuve grid search for parameters for TSNE and UMAP"""
import argparse
import Bio.SeqIO
from collections import OrderedDict
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import re
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
    parser.add_argument("--distance-matrix", help="csv file with the distance matrix")
    parser.add_argument("--alignment", help="fasta file with the alignment")
    parser.add_argument("--node-data", help="csv file with the clade_membership - that MUST be the name of the column.")
    parser.add_argument("--n-repeats", type=int, help="the number of times the k fold generator should repeat the data")
    parser.add_argument("--output", help="the csv path where the best thresholds will be saved.")
    parser.add_argument("--output-metadata", help="the path where the cross validation data will be saved.")
    parser.add_argument("--output-figure", help="PNG with the results displayed graphically")
    parser.add_argument("--output-total-data", help="output all the data to visualize training and visualization data")

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

    default_tuned_values = []
    list_of_embedding_strings = ["t-SNE","UMAP","MDS", "PCA"]
    embedding_class = [TSNE, UMAP, MDS, PCA]

    embedding_parameters = {
        "metric": "precomputed",
        "perplexity": 30,
        "learning_rate": 500.0,
        "square_distances" : True
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "init": "spectral",
        'n_neighbors': 200, 
        "min_dist": 0.05
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "dissimilarity": "precomputed",
        "n_components" : 2
    }
    default_tuned_values.append(embedding_parameters)

    embedding_parameters = {
        "n_components" : 10,
        "svd_solver" : "full"
    }
    default_tuned_values.append(embedding_parameters)

    # reading in the distance matrix and node data
 
    distance_matrix = pd.read_csv(args.distance_matrix, index_col=0)
    strain = distance_matrix.index.values.tolist()
    #distance_matrix.reset_index(drop=True)
    #distance_matrix = distance_matrix.to_numpy()

    
    node_df = pd.read_csv(args.node_data, sep="\t")
    clade_annotations = node_df[["strain", "clade_membership"]]
    strains_df = pd.DataFrame(strain, columns=["strain"])
    clade_annotations = clade_annotations.merge(strains_df, on="strain")
    node_df = node_df.merge(strains_df, on="strain")

    distance_matrix.columns = distance_matrix.index
    indices_to_drop = distance_matrix[~distance_matrix.index.isin(clade_annotations["strain"])].dropna(how = 'all')
    distance_matrix = distance_matrix[distance_matrix.index.isin(clade_annotations["strain"])].dropna(how = 'all')
    distance_matrix = distance_matrix.drop(indices_to_drop.index, axis=1)
    sequence_names = distance_matrix.index.values.tolist()
    distance_matrix = distance_matrix.to_numpy()

    #sequence_names = node_df["strain"].values.tolist()
    # getting PCA values

    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(args.alignment, "fasta"):
        if sequence.id in sequence_names:
            sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names_val = list(sequences_by_name.keys())
    #import pdb; pdb.set_trace()
    assert(len(sequence_names_val) == len(sequence_names))

    numbers = list(sequences_by_name.values())[:]
    for i in range(0,len(list(sequences_by_name.values()))):
        numbers[i] = re.sub(r'[^AGCT]', '5', numbers[i])
        numbers[i] = list(numbers[i].replace('A','1').replace('G','2').replace('C', '3').replace('T','4'))
        numbers[i] = [int(j) for j in numbers[i]]

    numbers = np.array(numbers)

    if args.output_total_data is not None:
        visualize_df = pd.DataFrame()

    #cross validation
    random_state = 12883823
    rkf = RepeatedKFold(n_splits=2,  n_repeats=args.n_repeats, random_state=random_state)
    k = 0
    total_list_methods = []

    for training_index, validation_index in rkf.split(sequence_names): 
        i = 0
        print("here " + str(k))
        for embed in default_tuned_values:
            print(i)
        # Calculate Euclidean distance between pairs of samples in the embedding.
        # The output is a condensed distance matrix with distances between pairs.
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
            scaler = StandardScaler()
            if(list_of_embedding_strings[i] == "PCA"):
                #performing PCA on my pandas dataframe
                numbers_subset = numbers[training_index]
                print(numbers_subset)
                print(set(training_index))
                pca = PCA(**embed) #can specify n, since with no prior knowledge, I use None
                training_embedding = pca.fit_transform(numbers_subset)
            else:
                # Subset distance matrix to training indices.
                training_distance_matrix = distance_matrix[training_index][:, training_index]

                # Embed training distance matrix.
                embedder = embedding_class[i](**embed)
                training_embedding = embedder.fit_transform(training_distance_matrix)
                

            # Calculate Euclidean distance between pairs of samples in the embedding.
            # The output is a condensed distance matrix with distances between pairs.
            # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
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
            genetic_classifier = LinearSVC(
                dual=False,
                random_state=0,
                class_weight={1: 5},
                verbose=0
            )
            X = scaler.fit_transform(squareform(training_distance_matrix).reshape(-1, 1))
            y = training_clade_status_for_pairs
            genetic_classifier.fit(X, y)

            x_range = np.linspace(-3, 3, 1000)
            z = genetic_classifier.decision_function(x_range.reshape(-1, 1))
            classifier_threshold_genetic = x_range[np.argwhere(z > 0)[-1]][0]
            
            if(list_of_embedding_strings[i] == "PCA"):
                #performing PCA on my pandas dataframe
                numbers_subset = numbers[validation_index]
                pca = PCA(**embed) #can specify n, since with no prior knowledge, I use None
                validation_embedding = pca.fit_transform(numbers_subset)
            else:
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

            if args.output_total_data is not None:
                if list_of_embedding_strings[i] != "PCA":
                    training_embedding_df = pd.DataFrame(training_embedding, columns=[str(list_of_embedding_strings[i])+"_x" + str(k), str(list_of_embedding_strings[i]) + "_y" + str(k)])
                    validation_embedding_df = pd.DataFrame(validation_embedding, columns=[str(list_of_embedding_strings[i])+"_x" + str(k), str(list_of_embedding_strings[i]) + "_y" + str(k)])
                else:
                    training_embedding_df = pd.DataFrame(training_embedding, columns=["PCA" + str(k) + str(j) for j in range(1,11)])
                    validation_embedding_df = pd.DataFrame(validation_embedding, columns=["PCA" + str(k) + str(j) for j in range(1,11)])
                
                training_embedding_df["clade_membership"] = clade_annotations.loc[training_index, "clade_membership"].values
                validation_embedding_df["clade_membership"] = clade_annotations.loc[validation_index, "clade_membership"].values

                merged_df = validation_embedding_df.merge(training_embedding_df, how='outer', left_index=True, right_index=True)
                visualize_df = visualize_df.append(merged_df)
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

            matthews_cc_val_genetic = matthews_corrcoef(validation_clade_status_for_pairs, genetic_classifier.predict(validation_embedding_distances))

            method_dict = {}
            
            method_dict["method"] = list_of_embedding_strings[i]

            method_dict["confusion_matrix"] = confusion_matrix_val
            method_dict["matthews_cc"] = matthews_cc_val
            method_dict["threshold"] = classifier_threshold
            method_dict["accuracy"] = accuracy


            print(method_dict)
            total_list_methods.append(method_dict)
            i = i + 1
        
        total_list_methods.append({"method":"genetic", "matthews_cc": matthews_cc_val_genetic, "threshold": classifier_threshold_genetic})
        print({"method":"genetic", "matthews_cc": matthews_cc_val_genetic, "threshold": classifier_threshold_genetic})
        k = k + 1


    cross_v_info = pd.DataFrame(total_list_methods)

    print(cross_v_info)
    if args.output_metadata is not None: 
        cross_v_info.to_csv(args.output_metadata)

    if args.output_figure is not None:

        cross_v_means = cross_v_info.groupby("method")["matthews_cc"].mean().reset_index()
        average = []
        for i in cross_v_info["method"].unique().tolist():
            cross_v_info_method = cross_v_info.where(cross_v_info["method"] == i).dropna(how = "all")
            average.append(np.average(np.array(cross_v_info_method["matthews_cc"].values.tolist()).flatten()))
        
        average.append(np.average(cross_v_info[cross_v_info['method'] == 'genetic']["matthews_cc"].values.tolist()))
        
        average = dict(zip(cross_v_info["method"].unique().tolist(), average))

        df = pd.DataFrame.from_dict(average, orient="index")
        df["method"] = df.index
        df.columns=['mean', 'method']
        df = df.reset_index(drop=True)
        #cross_v_info = cross_v_info.append(df)


        mpl.style.use("seaborn")
        plt.scatter(x=cross_v_info["method"], y=cross_v_info["matthews_cc"])
        plt.scatter(x=cross_v_means["method"], y=cross_v_means["matthews_cc"], marker="_")
        #plt.scatter(x=cross_v_info["method"], y=cross_v_info["mean"], marker="_")
        plt.ylim(0, 1)

        plt.savefig(args.output_figure)

    if args.output is not None:
        list_of_best_method = []
        list_of_best_threshold = []
        for i in ['PCA', 'MDS', 't-SNE', 'UMAP', 'genetic']:
            df = cross_v_info[cross_v_info.method == i].dropna(how="all")
            print(df)
            val = np.average(df["threshold"].values.tolist())
            print("val is below")
            print(val)
            list_of_best_method.append(i)
            list_of_best_threshold.append(val)
        output_df = pd.DataFrame(zip(list_of_best_method, list_of_best_threshold))
        output_df.columns = ["method", "threshold"]
        output_df.to_csv(args.output)

    if args.output_total_data is not None:
        visualize_df.to_csv(args.output_total_data)