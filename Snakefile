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

DISTANCE_THRESHOLDS = [
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0,
    5.5,
    6.0,
]
CLUSTER_MIN_SIZE = 10
CLUSTER_MIN_SAMPLES = 5

INTERNAL_NODE = [
    "sequences",
]

RANDOM_SEED = 314159
CLUSTER_THRESHOLD = 2.0

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
include: "sars-cov-2-nextstrain/Snakefile"
include: "sars-cov-2-nextstrain-2022-2023/Snakefile"

rule pathogens:
    input:
        *rules.seasonal_flu_training.input,
        *rules.seasonal_flu_test.input,
        *rules.seasonal_flu_reassortment.input,
        *rules.sarscov2.input,

# Include rules for the manuscript.
include: "docs/Snakefile"
