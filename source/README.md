# Cartography
Reduced dimension embeddings for pathogen sequences

Cartography is an open-source software for scientists, epidemiologists, etc. to run reduced dimension embeddings (PCA, MDS, t-SNE, and UMAP) on viral populations. This is the source code from the paper Cartography written by Sravani Nanduri and John Huddleston.

Source Code:
Bug reports:
Paper:


## Installing the package

Simply install the package using pip.

'''
pip install CartographyGen
'''
## Using Cartography

Cartography can be used via shell command:

'''
python3 CartographyGen \
    --distance-matrix {input.distance_matrix} \
    --alignment {input.alignment} \
    --cluster-data {input.cluster} \
    --cluster-threshold {input.threshold} \
    --random-seed {input.random-seed} \
    --output-node-data {output.node_data} \
    --output-dataframe {output.dataframe} \
    --output-figure {output.figure} \
    method=pca,mds,t-sne,umap \
    PCA:
        --components
        --explained-variance
    MDS:
        --components {params.components} \
    t-SNE:
        --perplexity
        --learning-rate
    UMAP:
    --nearest-neighbors
    --min-dist
'''

## Parameters:

- Distance Matrix: a csv distance matrix that can be read in by pandas, index column as row 0
- Alignment
- Cluster Data
- Cluster Threshold
- Random Seed
- Output-node-data
- Output-dataframe
- Output-figure

#### Method Level Parameters

PCA: 

