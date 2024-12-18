#!/usr/bin/env python3

import itertools
import pandas as pd

methods = ["pca", "mds", "t-sne", "umap"]
distance_thresholds = range(0, 16, 2)
n_components = [2,4,6]

# t-SNE
learning_rates = [100, 200, 500]
perplexity = [15, 30, 100]
columns = (
    "distance_threshold",
    "perplexity",
    "learning_rate",
)
tsne_parameters = itertools.product(
    distance_thresholds,
    perplexity,
    learning_rates
)
tsne_df = pd.DataFrame(
    tsne_parameters,
    columns=columns
)
tsne_df["method"] = "t-sne"

# UMAP
# Smaller values focus on local structure over global
n_neighbors = [25, 50, 100]
# Smaller values allow points to cluster closer together.
min_dist = [0.05, 0.1, 0.25]
columns = (
    "distance_threshold",
    "n_neighbors",
    "min_dist",
)
umap_parameters = itertools.product(
    distance_thresholds,
    n_neighbors,
    min_dist
)
umap_df = pd.DataFrame(
    umap_parameters,
    columns=columns
)
umap_df["method"] = "umap"

# MDS
columns = (
    "distance_threshold",
    "n_components"
)

mds_parameters = itertools.product(
     distance_thresholds,
     n_components
)
mds_df = pd.DataFrame(
    mds_parameters,
    columns=columns
)
mds_df["method"] = "mds"

# PCA
pca_df = pd.DataFrame({
    "method": "pca",
    "distance_threshold": distance_thresholds
})

# Collect methods parameters.
df = pd.concat([
    pca_df,
    mds_df,
    tsne_df,
    umap_df,
], ignore_index=True)

# Maintain the same column order.
columns = (
    "distance_threshold",
    "n_components",
    "perplexity",
    "learning_rate",
    "method",
    "n_neighbors",
    "min_dist",
)
df["n_components"] = df.n_components.astype('Int64')


df.to_csv(
    #"../sars-cov-2-nextstrain/config/method_parameters.tsv",
    sep="\t",
    index=False,
    na_rep="N/A",
    columns=columns,
)