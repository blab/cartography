from snakemake.utils import min_version
min_version("6.0")

module seasonal_flu_training_workflow:
    snakefile:
        "seasonal-flu-nextstrain/Snakefile"

use rule * from seasonal_flu_training_workflow as seasonal_flu_training_*
