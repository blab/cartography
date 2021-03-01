import argparse
import altair as alt
from altair_saver import save
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

#Figures to generate: MDS with LOF circles colored by FN, TP, etc ; 1D distance LOF scores with predicted and true outlier status; 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", help= "name of embedding type")
    parser.add_argument("--embedding", help="csv file with the embedding infomation")
    parser.add_argument("--true-outliers", help="the true outlier values. If given alone with find-outlier, the confusion matrix, matthews correlation coefficient, etc will be given")
    parser.add_argument("--find-outlier", action="store_true", help="no true outlier information is given, find outliers only")
    parser.add_argument("--columns", nargs=2, help="2 columns in the data to plot and find outlier in")
    parser.add_argument("--output-outliers", required=False, help="list of outliers")
    parser.add_argument("--output-main-figure", help="PNG with outlier circles describring the LOF score of the points in the graph")
    parser.add_argument("--output-LOF-figure", help="PNG with LOF distances colored by outlier and predicted outlier status")
    parser.add_argument("--output-metadata", help="the path where the outlier accuracy information should be saved. This will only work if both find-outlier and true-outliers are defined.")

    args = parser.parse_args()

    if args.output_metadata is not None and (args.true_outliers is None or args.find_outlier is None):
        print("you must provide true outliers and select find-outlier in order to output metadata about accuracy.", file=sys.stderr)
        sys.exit(1)

    if args.output_LOF_figure is not None and (args.find_outlier):
        print("you must provide true outliers for this figure.", file=sys.stderr)
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

    predicted = clf.fit_predict(total_list)
    X_scores = clf.negative_outlier_factor_

    if args.output_outliers is not None:
        
        classifier_threshold = (np.mean(X_scores) - (4*np.std(X_scores)))
        
        estimated_outlier_status = np.where(X_scores < classifier_threshold, -1, 1)

        embedding_df["predicted_outlier_status"] = estimated_outlier_status

        embedding_df["predicted_LOF_outlier_status"] = predicted

        embedding_df["X_scores"] = X_scores

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

    

    if args.output_main_figure is not None:

        plt.title("Local Outlier Factor (LOF)")
        if args.find_outlier:
            predicted_outliers = embedding_df["predicted_outlier_status"].values.tolist()
            confusion_matrix_values = []
            for i in range(len(predicted)): 
                #Not Outlier
                if predicted_outliers[i]==1:
                    confusion_matrix_values.append('#0000FF')
                #Outlier
                elif predicted_outliers[i]==-1:
                    confusion_matrix_values.append('#FF6600')
            from matplotlib.lines import Line2D
            legend_elements = [Line2D([0], [0], color='#0000FF', lw=4, label='Not Outlier'),
                            Line2D([0], [0], color='#FF6600', lw=4, label='Outlier')]
        else:
            from matplotlib.lines import Line2D
            legend_elements = [Line2D([0], [0], color='#0000FF', lw=4, label='True'),
                            Line2D([0], [0], color='#FF6600', lw=4, label='False')]
            #creating the true/false positives/negative values for this
            
            true_outliers = embedding_df["outlier"].values.tolist()
            predicted_outliers = embedding_df["predicted_outlier_status"].values.tolist()
            predicted = [int(x) for x in predicted_outliers]
            confusion_matrix_values = []
            for i in range(len(predicted)): 
                #TP
                if true_outliers[i]==1 and predicted[i]==1:
                    confusion_matrix_values.append('#0000FF')
                #FP
                elif predicted[i]==1 and true_outliers[i]== -1:
                    confusion_matrix_values.append('#FF6600')
                #TN
                elif true_outliers[i]==-1 and predicted[i]==-1:
                    confusion_matrix_values.append('#0000FF')
                #FN
                elif predicted[i]==-1 and true_outliers[i]==1:
                    confusion_matrix_values.append('#FF6600')
                else:
                    print(str(predicted[i]) + " " + str(true_outliers[i]))
            
        embedding_df["confusion_matrix_values"] = confusion_matrix_values
        plt.scatter(embedding_df[args.columns[0]].values.tolist(), embedding_df[args.columns[1]].values.tolist(), s=3., c=confusion_matrix_values)

        # plot circles with radius proportional to the outlier scores
        radius = (X_scores.max() - X_scores) / (X_scores.max() - X_scores.min())
        plt.scatter(embedding_df[args.columns[0]].values.tolist(), embedding_df[args.columns[1]].values.tolist(), s= 1000 * radius, edgecolors="#999999",
                    facecolors='none', label='Outlier scores')

        plt.axis('tight')
        plt.xlim((-250, 250))
        plt.ylim((-250, 300))

        legend = plt.legend(handles=legend_elements, loc="upper left")
        plt.savefig(args.output_main_figure)

    if args.output_LOF_figure is not None:
        embedding_df["LOF_scores"] = X_scores
        chart1 = alt.Chart(embedding_df).mark_circle(size=60).encode(
            x=alt.X('LOF_scores', title="Outliers"),
            color='outlier:N',
            tooltip=['strain']
        ).interactive()


        chart2 = alt.Chart(embedding_df).mark_circle(size=60).encode(
            x=alt.X('LOF_scores', title="Predicted Outliers"),
            color='predicted_outlier_status:N',
            tooltip=['strain', "outlier:N"]
        ).interactive()

        full = chart1|chart2
        full.save(args.output_LOF_figure)

