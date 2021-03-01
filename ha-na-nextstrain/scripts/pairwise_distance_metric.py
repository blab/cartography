"""performs procrustes analysis on the two embeddings given, calculates distance between them, returns values as a pandas dataframe. Can also return a procrustes analysis figure for you (if clade membership is given, it will be colored by that"""
import argparse
from augur.utils import read_node_data
from augur.utils import write_json
import numpy as np
import pandas as pd
from scipy.spatial import procrustes
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", help="name of embedding type")
    parser.add_argument("--embeddings", nargs=2, help="embeddings to perform procrustes analysis on")
    parser.add_argument("--columns", nargs=2, help="names of columns in both embeddings to use in the procrustes analysis")
    parser.add_argument("--metadata", help="the node data for clades membership for the figure")
    parser.add_argument("--colors", nargs="+", help="a list of the colors to use in the analysis")
    parser.add_argument("--domain", nargs="+", help="the clade names")
    parser.add_argument("--output-distance-metric", help="JSON file of the distance for exporting to the tree")
    parser.add_argument("--output-figure",  help="PNG figure of the procrustes analysis with lines connecting the points. If clade membership is given in the metadata, it will be used")
    parser.add_argument("--output-boxplot",  help="PNG figure of the distances for boxplot split by the first embeddings clade membership")
    parser.add_argument("--output-metadata",  help="extra information about the data given")

    args = parser.parse_args()

    if args.output_distance_metric is None and args.output_boxplot is not None:
        print("You must create the distance metric to create the boxplot", file=sys.stderr)
        sys.exit(1)

    if args.metadata is None and args.output_boxplot is not None:
        print("You must have metadata to create the boxplot", file=sys.stderr)
        sys.exit(1)

    embedding_1_df = pd.read_csv(args.embeddings[0])
    embedding_2_df = pd.read_csv(args.embeddings[1])

    if args.metadata is not None:
        node_data = read_node_data(args.metadata)
        metadata_df = clade_annotations = pd.DataFrame([
            {"strain": strain, "clade_membership": annotations["clade_membership"]}
            for strain, annotations in node_data["nodes"].items()])
        embedding_1_df = metadata_df.merge(embedding_1_df, on="strain")
        embedding_2_df = metadata_df.merge(embedding_2_df, on="strain")

    #procrustes analysis on the embeddings
    a = np.array([list(a) for a in zip(embedding_1_df[args.columns[0]].values.tolist(), embedding_1_df[args.columns[1]].values.tolist())])
    b = np.array([list(a) for a in zip(embedding_2_df[args.columns[0]].values.tolist(), embedding_2_df[args.columns[1]].values.tolist())])
    mtx1, mtx2, disparity = procrustes(a, b)

    df1 = pd.DataFrame(mtx1.tolist(), columns =["1_scaled_x", "1_scaled_y"])  
    df2 = pd.DataFrame(mtx2.tolist(), columns =["2_scaled_x", "2_scaled_y"])  
    merged_scaled_df = pd.merge(df1, df2, left_index=True, right_index=True)
    merged_scaled_df["strain"] = embedding_1_df["strain"]
    merged_scaled_df["clade_membership"] = embedding_1_df["clade_membership"]

    if args.output_distance_metric is not None:
        Ax = merged_scaled_df["1_scaled_x"].to_numpy()
        Ay = merged_scaled_df["1_scaled_y"].to_numpy()
        Ox = merged_scaled_df["2_scaled_x"].to_numpy()
        Oy = merged_scaled_df["2_scaled_y"].to_numpy()

        distance = np.sqrt(np.sum(((Ax-Ox)**2, (Ay-Oy)**2), axis=0))
        merged_scaled_df["distance"] = distance

        if args.output_metadata is not None:
            merged_scaled_df.to_csv(args.output_metadata)

        classifier_threshold = (np.mean(distance) + (1*np.std(distance)))
        estimated_outlier_status = np.where(distance < classifier_threshold, -1, 1)

        distance_df = pd.DataFrame()
        distance_df["distance_" + str(args.method)] = estimated_outlier_status

        #distance_df["distance_" + str(args.method)] = distance
        distance_df.index = merged_scaled_df["strain"]
        distance_dict = distance_df.transpose().to_dict()
        write_json({"nodes": distance_dict}, args.output_distance_metric)


        if args.output_boxplot is not None:
            sns_plot = sns.catplot(x="clade_membership", y="distance", kind="box", data=merged_scaled_df, height=4, aspect = 2)
            sns_plot.savefig(args.output_boxplot)

    if args.output_figure is not None:
        if args.metadata is not None:
            from matplotlib.lines import Line2D
            domain = args.domain
            range_ = args.colors
            print(range_)
            legend_elements = [Line2D([0], [0], color=range_[i], lw=4, label=domain[i]) for i in range(0, len(domain))]
            print(legend_elements)
            df = merged_scaled_df.copy()
            df.replace(dict(zip(domain,range_)), inplace=True)
            val = df["clade_membership"].to_numpy()
            print(val)
        
        line_segments = []
        for i in range(0, len(mtx1)):
            line_segments.append([(mtx1[i,0], mtx1[i,1]), (mtx2[i,0], mtx2[i,1])])

        x = merged_scaled_df["1_scaled_x"].values
        y = merged_scaled_df["1_scaled_y"].values

        pos_x = merged_scaled_df["2_scaled_x"].values
        pos_y = merged_scaled_df["2_scaled_y"].values

        fig, ax = plt.subplots(figsize=(8,8))
        if args.metadata is not None:
            ax.scatter(pos_x, pos_y, c=val)
            ax.legend(handles=legend_elements, loc=1)
        else:
            ax.scatter(pos_x, pos_y)

        ax.scatter(x, y, c="#696969")
        if args.metadata is not None:
            line_segments = LineCollection(line_segments, linewidths=(0.5, 1, 1.5, 2), color=val, linestyle='solid', alpha=0.5)
        else:
            line_segments = LineCollection(line_segments, linewidths=(0.5, 1, 1.5, 2), linestyle='solid', alpha=0.5)
        
        ax.add_collection(line_segments)
        ax.set_xlim(min(pos_x) - .015 , max(pos_x) + .015)
        ax.set_ylim(min(pos_y) - .015 , max(pos_y) + .015)

        plt.savefig(args.output_figure, dpi=300)