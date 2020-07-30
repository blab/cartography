"""
takes tree path, disease name, and similarity matrix and outputs node_df and re-done similarity_matrix
"""
import pandas as pd
import numpy as np
from augur.utils import json_to_tree
import json
from pathlib import Path
from Helpers import get_y_positions
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--tree_path", required=True, help="tree to convert for altair purposes")
parser.add_argument("--disease_name", required=True, help="name of disease")
	
args = parser.parse_args()

similarity_matrix = pd.read_csv('../Dataframes/distance_matrix_' + args.disease_name)

with open(args.tree_path) as fh:
    json_tree_handle = json.load(fh)

tree = json_to_tree(json_tree_handle)

heights = get_y_positions(tree)
for node in tree.find_clades():
    node.yvalue = heights[node]
	
node_data = [
    {
        "strain": node.name,
        "date": node.node_attrs["num_date"]["value"],
        "y": node.yvalue,
        "region": node.node_attrs["region"]["value"],
        "country": node.node_attrs["country"]["value"],
        "parent_date": node.parent is not None and node.parent.node_attrs["num_date"]["value"] or node.node_attrs["num_date"],
        "parent_y": node.parent is not None and node.parent.yvalue or node.yvalue,
        "clade_membership" : node.node_attrs['clade_membership']["value"]
    }
    for node in tree.find_clades(terminal=True)
]

node_df = pd.DataFrame(node_data)

node_df["y"] = node_df["y"].max() - node_df["y"]

node_df["parent_y"] = node_df["parent_y"].max() - node_df["parent_y"]

# Reannotate clades that we aren't interested in as "other" to simplify color assignment in visualizations.
try:
    node_df["clade_membership_color"] = node_df["clade_membership"].apply(lambda clade: clade if clade in clades_to_plot else "other")
except:
    node_df["clade_membership_color"] = node_df["clade_membership"]
	
indices_to_drop = similarity_matrix[~similarity_matrix.index.isin(node_df["strain"])].dropna(how = 'all')
similarity_matrix = similarity_matrix[similarity_matrix.index.isin(node_df["strain"])].dropna(how = 'all')
similarity_matrix = similarity_matrix.drop(indices_to_drop.index, axis=1)

similarity_matrix.to_csv('../Dataframes/distance_matrix_' + args.disease_name + ".csv")

node_df.to_csv('../Dataframes/node_dataframe' + args.disease_name + '.csv')
	