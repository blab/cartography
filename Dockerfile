FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="46e372c7772af8a6dc79f1d5ef885e12edf75fe9839d9b38ee010668be95888c"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: cartography.yml
#   prefix: /conda-envs/4f5bdb38739416e451545e49f72c0b3d
#   name: cartography
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - augur=13.1.2
#     - pip
#     - seaborn
#     - altair
#     - altair_saver
#     - dendropy
#     - python=3.7*
#     - scikit-learn
#     - umap-learn
#     - jsonschema
#     - jupyterlab
#     - matplotlib
#     - nodejs
#     - pandoc=2.14.1
#     - pandoc-crossref=0.3.12.0
#     - selenium
#     - python-chromedriver-binary
#     - reportlab
#     - wget
#     - snakemake
#     - snp-sites
#     - tabix
#     - pixy
#     - openjdk=11
#     - tsv-utils
#     - samtools
#     - joblib=1.1.0 # required for hdbscan to work
#     - hdbscan=0.8.28
#     - statsmodels
#     - pip:
#       - pathogen-embed==0.1.0
#       - statistics
#       - svglib
#       - tabulate
RUN mkdir -p /conda-envs/4f5bdb38739416e451545e49f72c0b3d
COPY cartography.yml /conda-envs/4f5bdb38739416e451545e49f72c0b3d/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/4f5bdb38739416e451545e49f72c0b3d --file /conda-envs/4f5bdb38739416e451545e49f72c0b3d/environment.yaml && \
    mamba clean --all -y

RUN apt-get update && \
    apt-get install -y wget && \
    rm -rf /var/lib/apt/list/*

ENV JULIA_VERSION=1.8.5

RUN mkdir /opt/julia-${JULIA_VERSION} && \
    cd /tmp && \
    wget -q https://julialang-s3.julialang.org/bin/linux/x64/`echo ${JULIA_VERSION} | cut -d. -f 1,2`/julia-${JULIA_VERSION}-linux-x86_64.tar.gz && \
    tar xzf julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C /opt/julia-${JULIA_VERSION} --strip-components=1 && \
    rm /tmp/julia-${JULIA_VERSION}-linux-x86_64.tar.gz

RUN ln -fs /opt/julia-*/bin/julia /usr/local/bin/julia
RUN apt-get update && apt-get install -y gcc g++ && rm -rf /var/lib/apt/lists/*
RUN julia -e 'using Pkg; Pkg.add(["TreeKnit"]);'

