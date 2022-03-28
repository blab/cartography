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
    "sequences",
]

RANDOM_SEED = 314159

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

# Include common functions and rules shared across analyses.
include: "rules/common.smk"

# Include rules for each pathogen.
include: "seasonal-flu-nextstrain/Snakefile"
include: "seasonal-flu-nextstrain-2018-2020/Snakefile"
include: "ha-na-nextstrain/Snakefile"
include: "mers-nextstrain/Snakefile"
include: "mers-muller/Snakefile"
include: "sars-cov-2-nextstrain/Snakefile"

# include: "outlier_analysis/Snakefile"
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
#include: "docs/Snakefile"
