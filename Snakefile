from snakemake.utils import min_version, Paramspace
min_version("6.0")

import pandas as pd

# Define top-level configuration parameters.
EMBEDDING_METHODS = [
    "pca",
    "mds",
    "t-sne",
    "umap"
]

INTERNAL_NODE = [
    "ancestral",
    "sequences"
]

RANDOM_SEED = 314159

wildcard_constraints:
    method="(pca|mds|t-sne|umap)",
    internal_node= "(ancestral|sequences)"

# Define final outputs for the workflow.
rule all:
    input:
        # Static version of the paper.
        "docs/cartography.pdf",
        # Interactive version of the paper.
        "docs/cartography.html",

# Include rules for each pathogen.
include: "seasonal-flu-nextstrain/Snakefile"
# include: "seasonal-flu-nextstrain-2018-2020/Snakefile"
# include: "ha-na-nextstrain/Snakefile"
# include: "ha-na-ma-nextstrain/Snakefile"
# include: "outlier_analysis/Snakefile"
# include: "zika-nextstrain/Snakefile"
# include: "mers-nextstrain/Snakefile"
# include: "sars-cov-2-nextstrain/Snakefile"

# Include rules for the manuscript.
#include: "docs/Snakefile"
