#!/usr/bin/env python3
from augur.utils import read_node_data
import hdbscan
import matplotlib.pyplot as plt
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


def get_pairwise_clade_status(clades):
    """Traverse pairs of samples from left-to-right, top-to-bottom along the upper
    triangle of the pairwise matrix and collect the clade status of each pair as
    either within- or between-clades. This traversal excludes self-self
    comparisons along the diagonal.

    """
    clade_status = []
    for i in range(len(clades) - 1):
        for j in range(i + 1, len(clades)):
            if clades[i] == clades[j]:
                clade_status.append(1)
            else:
                clade_status.append(0)

    return np.array(clade_status)


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
    for key, value in snakemake.params.method_parameters.items()
    if value != "NA"
}
method_parameters.update(DEFAULT_PARAMETERS_BY_METHOD[method])
print(method_parameters)

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

    # Name columns to match the strains in the index.
    input_matrix.columns = input_matrix.index.values

    # Reorder its rows to match the clades above.
    input_matrix = input_matrix.loc[strains, strains]
    input_matrix_strains = input_matrix.index.values
    assert np.array_equal(strains, input_matrix_strains)

    # Extract the numpy arrays corresponding to the distance matrix data frame.
    input_matrix = input_matrix.values
    is_distance_matrix = True

print(input_matrix)

folds = RepeatedKFold(
    n_splits=N_SPLITS,
    n_repeats=N_REPEATS,
    random_state=RANDOM_STATE
)

all_results = []
for cv_iteration, (training_index, validation_index) in enumerate(folds.split(strains)):
    print(f"Iteration: {cv_iteration}")
    results = method_parameters.copy()
    results["method"] = method
    results["cv_iteration"] = cv_iteration
    results["distance_threshold"] = distance_threshold

    training_matrix = input_matrix[training_index]
    training_clades = clades.iloc[training_index]["clade_membership"].values
    training_observed_clade_status = get_pairwise_clade_status(training_clades)

    validation_matrix = input_matrix[validation_index]
    validation_clades = clades.iloc[validation_index]["clade_membership"].values
    validation_observed_clade_status = get_pairwise_clade_status(validation_clades)

    # Distance matrix is square and needs to be sliced two ways.
    if is_distance_matrix:
        training_matrix = training_matrix[:, training_index]
        validation_matrix = validation_matrix[:, validation_index]

    # Train the embedding method.
    embedder = method_class(**method_parameters)
    training_embedding = embedder.fit_transform(training_matrix)

    plt.plot(training_embedding[:, 0], training_embedding[:, 1], "o")
    plt.savefig(snakemake.output.table.replace(".tsv", f"_{cv_iteration}_training_embedding.pdf"))
    plt.close()

    # Calculate fit of clustering to trained embedding.
    clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=distance_threshold)
    clusterer.fit(training_embedding)
    clusters = clusterer.labels_.astype(str)
    training_predicted_clade_status = get_pairwise_clade_status(clusters)

    training_mcc = matthews_corrcoef(
        training_observed_clade_status,
        training_predicted_clade_status
    )
    results["training_mcc"] = training_mcc

    training_confusion_matrix = confusion_matrix(
        training_observed_clade_status,
        training_predicted_clade_status
    )
    results["training_tn"] = training_confusion_matrix[0, 0]
    results["training_fn"] = training_confusion_matrix[1, 0]
    results["training_tp"] = training_confusion_matrix[1, 1]
    results["training_fp"] = training_confusion_matrix[0, 1]

    # Validate the embedding method.
    embedder = method_class(**method_parameters)
    validation_embedding = embedder.fit_transform(validation_matrix)

    plt.plot(validation_embedding[:, 0], validation_embedding[:, 1], "o")
    plt.savefig(snakemake.output.table.replace(".tsv", f"_{cv_iteration}_validation_embedding.pdf"))
    plt.close()

    # Calculate fit of clustering to trained embedding.
    clusterer = hdbscan.HDBSCAN(cluster_selection_epsilon=distance_threshold)
    clusterer.fit(validation_embedding)
    clusters = clusterer.labels_.astype(str)
    validation_predicted_clade_status = get_pairwise_clade_status(clusters)

    validation_confusion_matrix = confusion_matrix(
        validation_observed_clade_status,
        validation_predicted_clade_status
    )
    results["validation_tn"] = validation_confusion_matrix[0, 0]
    results["validation_fn"] = validation_confusion_matrix[1, 0]
    results["validation_tp"] = validation_confusion_matrix[1, 1]
    results["validation_fp"] = validation_confusion_matrix[0, 1]

    results["validation_mcc"] = matthews_corrcoef(
        validation_observed_clade_status,
        validation_predicted_clade_status
    )

    all_results.append(results)

df = pd.DataFrame(all_results)
df.to_csv(
    snakemake.output.table,
    sep="\t",
    index=False
)
