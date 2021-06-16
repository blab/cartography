#!/usr/bin/env python3
from augur.utils import read_node_data
import hdbscan
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from sklearn.model_selection import RepeatedKFold
from umap import UMAP

import sys
sys.path.append("../notebooks/scripts/")

from Helpers import get_PCA_feature_matrix

N_SPLITS = 2
N_REPEATS = 3
RANDOM_STATE = 12883823
CLASS_BY_METHOD = {
    "pca": PCA,
    "mds": MDS,
    "t-sne": TSNE,
    "umap": UMAP,
}
DEFAULT_PARAMETERS_BY_METHOD = {
    "pca": {
        "n_components": 10,
        "svd_solver": "full",
    },
    "mds": {
        "dissimilarity": "precomputed",
        "n_components": 2,
        "n_init": 2,
        "n_jobs": 1,
    },
    "t-sne": {
        "metric": "precomputed",
        "square_distances": True,
    },
    "umap": {
        "init": "spectral",
    },
}


print(snakemake.input)
print(snakemake.output)
print(snakemake.params)

method = snakemake.params.method_parameters.pop("method")
method_class = CLASS_BY_METHOD[method]
distance_threshold = float(snakemake.params.method_parameters.pop("distance_threshold"))
method_parameters = {
    key: value
    for key, value in snakemake.params.method_parameters
    if value != "NA"
}
method_parameters.update(DEFAULT_PARAMETERS_BY_METHOD[method])

# Load clade annotations.
clades = read_node_data(snakemake.input.clades)["nodes"]
clades = [
    {"strain": strain, "clade_membership": values["clade_membership"]}
    for strain, values in clades.items()
    if not strain.startswith("NODE")
]
clades = pd.DataFrame(clades)
strains = clades["strain"].values
print(clades)
print(strains)

if method == "pca":
    # Load alignment.
    input_matrix = get_PCA_feature_matrix(snakemake.input.alignment, strains)
    is_distance_matrix = False
else:
    # Load distance matrix.
    input_matrix = pd.read_csv(snakemake.input.distance_matrix, index_col=0)

    # Reorder its rows to match the clades above.
    input_matrix = input_matrix.loc[strains]
    input_matrix_strains = input_matrix.index.values
    assert np.array_equal(strains, input_matrix_strains)

    # Extract the numpy arrays corresponding to the distance matrix data frame.
    input_matrix = input_matrix.values
    is_distance_matrix = True

print(input_matrix)

df = pd.DataFrame(snakemake.params)

folds = RepeatedKFold(
    n_splits=N_SPLITS,
    n_repeats=N_REPEATS,
    random_state=RANDOM_STATE
)

for cv_iteration, (training_index, validation_index) in enumerate(folds.split(strains)):
    print(f"Iteration: {cv_iteration}")
    training_matrix = input_matrix[training_index]
    training_clades = clades.iloc[training_index]["clade_membership"].values
    validation_matrix = input_matrix[validation_index]
    validation_clades = clades.iloc[validation_index]["clade_membership"].values

    # Distance matrix is square and needs to be sliced two ways.
    if is_distance_matrix:
        training_matrix = training_matrix[:, training_index]
        validation_matrix = validation_matrix[:, validation_index]

    # Train the embedding method.
    embedder = method_class(**method_parameters)
    training_embedding = embedder.fit_transform(training_matrix)

    # Calculate fit of clustering to trained embedding.
    clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=distance_threshold)
    clusterer.fit(training_embedding)
    clusters = clusterer.labels_.astype(str)

    #import ipdb; ipdb.set_trace()
    confusion_matrix_training = confusion_matrix(
        training_clades,
        clusters
    )
    matthews_correlation_coefficient = matthews_corrcoef(
        training_clades,
        clusters
    )
    print(f"confusion matrix: {confusion_matrix_training}")
    print(f"MCC: {matthews_correlation_coefficient}")

df.to_csv(
    snakemake.output.table,
    sep="\t",
    index=False
)
