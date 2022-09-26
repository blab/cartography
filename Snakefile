from snakemake.utils import min_version, Paramspace
min_version("6.0")

import pandas as pd

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

# Define top-level configuration parameters.
EMBEDDING_METHODS = [
    "pca",
    "mds",
    "t-sne",
    "umap"
]

INTERNAL_NODE = [
    "sequences",
]

RANDOM_SEED = 314159
CLUSTER_THRESHOLD = 2.0

localrules:
    seasonal_flu_training_aggregate_clusters_by_parameters,
    mers_download_elife_tree,
    mers_download_mcc_tree,
    mers_unzip,
    mers_aggregate_clusters_by_parameters,
    sarscov2_aggregate_clusters_by_parameters,

wildcard_constraints:
    method="(pca|mds|t-sne|umap|genetic)",
    internal_node= "(ancestral|sequences)",
    segment="(ha|na)",
    ha_concatenated="(ha|na|concatenated)",
    ha_concat="(ha|concatenated)"

# Define final outputs for the workflow.
rule all:
    input:
        # Static version of the paper.
        "docs/cartography.pdf",
        # Interactive version of the paper.
        "docs/cartography.html",

# Include rules for each pathogen.
include: "simulations/Snakefile"
include: "seasonal-flu-nextstrain/Snakefile"
include: "seasonal-flu-nextstrain-2018-2020/Snakefile"
include: "ha-na-nextstrain/Snakefile"
include: "mers-nextstrain/Snakefile"
include: "mers-muller/Snakefile"
include: "sars-cov-2-nextstrain/Snakefile"

# include: "zika-nextstrain/Snakefile"

rule pathogens:
    input:
        *rules.seasonal_flu_training.input,
        *rules.seasonal_flu_test.input,
        *rules.seasonal_flu_reassortment.input,
        *rules.mers.input,
        *rules.mers_muller.input,
        *rules.sarscov2.input,

# Include rules for the manuscript.
include: "docs/Snakefile"
