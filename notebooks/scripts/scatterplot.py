import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress
import seaborn as sns
import statsmodels.api
import statistics
import sys

from Helpers import scatterplot_xyvalues

if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description = "creates embeddings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--distance", required=True, help="a distance matrix of pairwise distances with the strain name as the index")
    parser.add_argument("--embedding", required=True, help="an embedding csv matrix - the order of distances per strain MUST be the same as the distance matrix")
    parser.add_argument("--method", required=True, choices = ["pca", "mds", "t-sne", "umap"], help="the embedding used")
    parser.add_argument("--output-figure", help="path for outputting as a PNG")
    parser.add_argument("--output-dataframe", help="path for outputting as a dataframe")
    
    args = parser.parse_args()
    
    #error handling
    
    if args.output_figure is None and args.output_dataframe is None:
        print("You must specify one of the outputs", file=sys.stderr)
        sys.exit(1)
        
    # reading in the distance matrix and embedding csv files, checking to make sure the format is correct
    
    distance_matrix = pd.read_csv(args.distance, index_col=0)
    embedding_df = pd.read_csv(args.embedding, index_col=0)
    assert np.array_equal(distance_matrix.index, embedding_df.index)

    #calling Helpers.py scatterplot_xyvalues on the data

    column_names = list(embedding_df.columns.values)

    total_df = scatterplot_xyvalues(list(embedding_df.index), distance_matrix, embedding_df, column_names[0], column_names[1], args.method)
        
    if args.output_figure is not None:
        y_values = statsmodels.nonparametric.smoothers_lowess.lowess(
        total_df["euclidean"],
        total_df["genetic"],
        frac=0.6666666666666666,
        it=3,
        delta=0.0,
        is_sorted=False,
        missing='drop',
        return_sorted=True
        )

        PD_Y_values = pd.DataFrame(y_values)
        PD_Y_values.columns = ["LOWESS_x", "LOWESS_y"]
        regression = linregress(total_df["genetic"], total_df["euclidean"])
        slope, intercept, r_value, p_value, std_err = regression
        
        if args.output_figure is not None:
            fig, ax = plt.subplots(1, 1, figsize=(6, 6))

            ax.plot(total_df["genetic"], total_df["euclidean"], "o", alpha=0.25)
            ax.plot(PD_Y_values["LOWESS_x"], PD_Y_values["LOWESS_y"], label="LOESS")

            ax.set_xlabel("Genetic distance")
            ax.set_ylabel(f"Euclidean distance ({args.method})")
            ax.set_title(f"Euclidean distance ({args.method}) vs. genetic distance ($R^2={r_value ** 2:.3f}$)")

            sns.despine()
            
            plt.savefig(args.output_figure)
            
    if args.output_dataframe is not None:
        total_df = pd.concat([total_df, PD_Y_values], axis=1)
        total_df["pearson_coef"]=r_value ** 2
        total_df.to_csv(args.output_dataframe)