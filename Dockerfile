# Start from a base image
FROM nextstrain/base:build-20240710T214955Z

# Install Python packages.
# Allow Snakemake to create subdirs in the user cache dir
# <https://github.com/nextstrain/ncov-ingest/pull/401>
RUN pip3 install \
    "altair[all]" \
    jupyterlab \
    notebook \
    seaborn \
    statsmodels \
 && rm -rf ~/.cache

# Install Java.
RUN apt-get update && apt-get install -y --no-install-recommends \
    default-jre && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Julia.

# Install TreeKnit.
