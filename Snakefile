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

EMBEDDING_NAME_BY_METHOD = {
    "pca": "PCA",
    "mds": "MDS",
    "t-sne": "t-SNE",
    "umap": "UMAP",
    "genetic": "genetic",
}

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
    6.5,
    7.0,
    7.5,
    8.0,
    8.5,
    9.0,
    9.5,
    10.0,
    10.5,
    11.0,
    11.5,
    12.0,
    12.5,
    13.0,
    13.5,
    14.0,
    14.5,
    15.0,
    15.5,
    16.0,
    16.5,
    17.0,
    17.5,
    18.0,
    18.5,
    19.0,
    19.5,
    20.0,
]
CLUSTER_MIN_SIZE = 10
CLUSTER_MIN_SAMPLES = 5

INTERNAL_NODE = [
    "sequences",
]

RANDOM_SEED = 314159

# Define parameters for replication of cluster accuracy analysis across late
# pathogen datasets.
CLUSTER_REPLICATION_REPLICATES = list(range(5))
CLUSTER_REPLICATION_SUBSAMPLE_MAX_SEQUENCES = (500, 1000, 1500, 2000, 2500)
CLUSTER_REPLICATION_SUBSAMPLING_SCHEMES = ("even", "random")

SEASONAL_FLU_REFERENCE_STRAIN = "A/Beijing/32/1992"
SARS_COV_2_REFERENCE_STRAIN = "Wuhan-Hu-1/2019"

wildcard_constraints:
    method="(pca|mds|t-sne|umap|genetic)",
    internal_node= "(ancestral|sequences)",
    segment="(ha|na)",
    ha_concatenated="(ha|na|concatenated)",
    ha_concat="(ha|concatenated)",
    clade_membership="(Nextstrain_clade|Nextclade_pango|Nextclade_pango_collapsed)",
    replicate="[0-9]+",
    subsampling_scheme="(even|random)",

# Define final outputs for the workflow.
rule all:
    input:
        # Static version of the paper.
        "manuscript/cartography.pdf",
        # Interactive version of the paper.
        #"manuscript/cartography.html",

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
        *rules.sarscov2_test.input,

# Include rules for the manuscript.
include: "manuscript/Snakefile"
