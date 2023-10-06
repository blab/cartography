# FROM condaforge/mambaforge:latest AS conda
# LABEL io.github.snakemake.containerized="true"
# LABEL io.github.snakemake.conda_env_hash="46e372c7772af8a6dc79f1d5ef885e12edf75fe9839d9b38ee010668be95888c"

# # Step 1: Retrieve conda environments

# # Conda environment:
# #   source: cartography.yml
# #   prefix: /conda-envs/4f5bdb38739416e451545e49f72c0b3d
# #   name: cartography
# #   channels:
# #     - conda-forge
# #     - bioconda
# #     - defaults
# #   dependencies:
# #     - augur=13.1.2
# #     - pip
# #     - seaborn
# #     - altair
# #     - altair_saver
# #     - dendropy
# #     - python=3.7*
# #     - scikit-learn
# #     - umap-learn
# #     - jsonschema
# #     - jupyterlab
# #     - matplotlib
# #     - nodejs
# #     - pandoc=2.14.1
# #     - pandoc-crossref=0.3.12.0
# #     - selenium
# #     - python-chromedriver-binary
# #     - reportlab
# #     - wget
# #     - snakemake
# #     - snp-sites
# #     - tabix
# #     - pixy
# #     - openjdk=11
# #     - tsv-utils
# #     - samtools
# #     - joblib=1.1.0 # required for hdbscan to work
# #     - hdbscan=0.8.28
# #     - statsmodels
# #     - pip:
# #       - pathogen-embed==0.1.0
# #       - statistics
# #       - svglib
# #       - tabulate
# RUN mkdir -p /conda-envs/4f5bdb38739416e451545e49f72c0b3d
# COPY cartography.yml /conda-envs/4f5bdb38739416e451545e49f72c0b3d/environment.yaml

# # Step 2: Generate conda environments

# RUN mamba env create --prefix /conda-envs/4f5bdb38739416e451545e49f72c0b3d --file /conda-envs/4f5bdb38739416e451545e49f72c0b3d/environment.yaml && \
#     mamba clean --all -y

FROM julia

RUN apt-get update -y && \
    apt-get install clang -y

RUN julia -e 'using Pkg; Pkg.add(["TreeKnit"]); Pkg.build("TreeKnit");'

ARG MINIFORGE_NAME=Miniforge3
ARG MINIFORGE_VERSION=23.1.0-4
ARG TARGETPLATFORM

ENV CONDA_DIR=/opt/conda
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_DIR}/bin:${PATH}

# 1. Install just enough for conda to work
# 2. Keep $HOME clean (no .wget-hsts file), since HSTS isn't useful in this context
# 3. Install miniforge from GitHub releases
# 4. Apply some cleanup tips from https://jcrist.github.io/conda-docker-tips.html
#    Particularly, we remove pyc and a files. The default install has no js, we can skip that
# 5. Activate base by default when running as any *non-root* user as well
#    Good security practice requires running most workloads as non-root
#    This makes sure any non-root users created also have base activated
#    for their interactive shells.
# 6. Activate base by default when running as root as well
#    The root user is already created, so won't pick up changes to /etc/skel
RUN apt-get update > /dev/null && \
    apt-get install --no-install-recommends --yes \
        wget bzip2 ca-certificates \
        git \
        tini \
        > /dev/null && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-$(uname -m).sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
    find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes  && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]

# Step 2: Generate conda environments

RUN conda upgrade conda

RUN conda install -c conda-forge mamba

COPY cartography.yml /cartography.yaml

RUN mamba env update -f cartography.yaml && \
    mamba clean --all -y
    
# RUN mamba env create --file cartography.yaml && \
#     mamba clean --all -y

#RUN conda activate cartography