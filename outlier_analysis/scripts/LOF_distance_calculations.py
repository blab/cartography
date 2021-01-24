import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from sklearn.neighbors import LocalOutlierFactor
from scipy.spatial import distance
from sklearn.svm import LinearSVC
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix, matthews_corrcoef, roc_curve
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", help= "name of embedding type")
    parser.add_argument("--embedding", help="csv file with the embedding infomation")
    parser.add_argument("--true-outliers", help="the true outlier values. If given alone with find-outlier, the confusion matrix, matthews correlation coefficient, etc will be given")
    parser.add_argument("--find-outlier", action="store_true", help="no true outlier information is given, find outliers only")
    parser.add_argument("--columns", nargs=2, help="2 columns in the data to plot and find outlier in")
    parser.add_argument("--output-outliers", required=False, help="list of outliers")
    parser.add_argument("--output-figure", help="PNG with outlier circles describring the LOF score of the points in the graph")
    parser.add_argument("--output-metadata", help="the path where the outlier accuracy information should be saved. This will only work if both find-outlier and true-outliers are defined.")

    args = parser.parse_args()

    if args.output_metadata is not None and (args.true_outliers is None or args.find_outlier is None):
        print("you must provide true outliers and select find-outlier in order to output metadata about accuracy.", file=sys.stderr)
        sys.exit(1)
    
    embedding_df = pd.read_csv(args.embedding)  
    
    if args.true_outliers is not None:
        metadata_df = pd.read_csv(args.true_outliers, sep="\t")

        embedding_df = embedding_df.merge(metadata_df, on="strain")

    mapping = {'not_outlier': 1, 'outlier': -1}

    embedding_df = embedding_df.replace({'not_outlier': mapping, 'outlier': mapping})

    clf = LocalOutlierFactor(n_neighbors=20, contamination='auto')

    print(args.columns[0])
    total_list = np.array([list(a) for a in zip(embedding_df['mds1'].values.tolist(), embedding_df['mds2'].values.tolist())])

    y_pred = clf.fit_predict(total_list)
    X_scores = clf.negative_outlier_factor_

    if args.output_outliers is not None:
        #create a SVM to find a distance threshold between the scores
        
        #initializing scaler
        scaler = StandardScaler()

        scaled_x_scores = scaler.fit_transform(X_scores.reshape(-1, 1)).flatten()

        print(scaled_x_scores)

        # Use a support vector machine classifier to identify an optimal threshold

        classifier = LinearSVC(
            dual=False,
            random_state=0,
            class_weight={-1: 5},
            verbose=0
        )
        X = np.array(scaled_x_scores).reshape(-1,1)
        y = embedding_df["outlier"].astype(float).values
        
        classifier.fit(X, y)

        # Find the SVM's threshold between the two given classes by passing the
        # range of possible scaled distance values (effectively z-scores) to the
        # classifier's decision function. See the documentation for more:
        # https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html#sklearn.svm.LinearSVC.decision_function
        
        x_range = np.linspace(-3, 3, 1000)
        z = classifier.decision_function(x_range.reshape(-1, 1))
        try:
            classifier_threshold = x_range[np.argwhere(z > 0)[-1]][0]
        except:
            classifier_threshold = "NaN"

        # Estimate group labels using the same input data used to train the SVM.
        estimated_outlier_status = classifier.predict(X)

        embedding_df["predicted_outlier_status"] = estimated_outlier_status

        if args.output_metadata is not None:
            values_df = pd.DataFrame()
            accuracy = accuracy_score(embedding_df["outlier"].values.tolist(), embedding_df["predicted_outlier_status"].values.tolist())
            matthews_cc_val = matthews_corrcoef(embedding_df["outlier"].values.tolist(), embedding_df["predicted_outlier_status"].values.tolist())
            confusion_matrix_val = confusion_matrix(embedding_df["outlier"].values.tolist(), embedding_df["predicted_outlier_status"].values.tolist())
            fpr, tpr, thresholds = roc_curve(embedding_df["outlier"].values.tolist(), embedding_df["predicted_outlier_status"].values.tolist())
            fpr_roc = fpr
            tpr_roc = tpr
            thresholds_roc = thresholds

            values_df = pd.DataFrame([[matthews_cc_val, accuracy, classifier_threshold, args.method, confusion_matrix_val[0][0], confusion_matrix_val[1][0], confusion_matrix_val[1][1], confusion_matrix_val[0][1], fpr_roc, tpr_roc, thresholds_roc]], columns=["MCC", "accuracy", "threshold", "embedding", "TN", "FN", "TP", "FP", "roc_fpr", "roc_tpr", "roc_thresholds"]).round(3)
            values_df.to_csv(args.output_metadata)
        

        embedding_df.to_csv(args.output_outliers, index=False)

    

    if args.output_figure is not None:

        plt.title("Local Outlier Factor (LOF)")
        if args.find_outlier:
            plt.scatter(embedding_df[args.columns[0]].values.tolist(), embedding_df[args.columns[1]].values.tolist(), color='k', s=3., label='Data points')
        else:
            groups = embedding_df.groupby('outlier')
            for name, group in groups:
                plt.scatter(embedding_df[args.columns[0]].values.tolist(), embedding_df[args.columns[1]].values.tolist(), s=3., label=name)

        # plot circles with radius proportional to the outlier scores
        radius = (X_scores.max() - X_scores) / (X_scores.max() - X_scores.min())
        plt.scatter(embedding_df[args.columns[0]].values.tolist(), embedding_df[args.columns[1]].values.tolist(), s= 1000 * radius, edgecolors='r',
                    facecolors='none', label='Outlier scores')

        plt.axis('tight')
        plt.xlim((-250, 250))
        plt.ylim((-250, 300))

        legend = plt.legend(loc='upper left')
        legend.legendHandles[0]._sizes = [15]
        legend.legendHandles[1]._sizes = [25]
        plt.savefig(args.output_figure)