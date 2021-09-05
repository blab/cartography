# Cartography
Reduced dimension embeddings for pathogen sequences


Cartography is an open-source software for scientists, epidemiologists, etc. to run reduced dimension embeddings (PCA, MDS, t-SNE, and UMAP) on viral populations. This is the source code from the paper Cartography written by Sravani Nanduri and John Huddleston.

Source Code:
Bug reports:
Paper:


## Installing the package

Simply install the package using pip.

```
pip install CartographyGen
```

# src.embed module

## Command line interface

Reduced dimension embeddings for pathogen sequences


```
usage: embed [-h] [--distance-matrix DISTANCE_MATRIX] [--separator SEPARATOR]
             [--alignment ALIGNMENT] [--cluster-data CLUSTER_DATA]
             [--cluster-threshold CLUSTER_THRESHOLD]
             [--random-seed RANDOM_SEED] [--output-dataframe OUTPUT_DATAFRAME]
             [--output-figure OUTPUT_FIGURE]
             {pca,t-sne,umap,mds} ...
```

### Positional Arguments

### Named Arguments

### Sub-commands:

#### pca

Principal Component Analysis

```
embed pca [-h] [--components COMPONENTS]
          [--explained-variance EXPLAINED_VARIANCE]
```

##### Named Arguments

#### t-sne

t-distributed Stochastic Neighborhood Embedding

```
embed t-sne [-h] [--perplexity PERPLEXITY] [--learning-rate LEARNING_RATE]
```

##### Named Arguments

#### umap

Uniform Manifold Approximation and Projection

```
embed umap [-h] [--nearest-neighbors NEAREST_NEIGHBORS] [--min-dist MIN_DIST]
```

##### Named Arguments

#### mds

Multidimensional Scaling

```
embed mds [-h] [--components COMPONENTS]
```

##### Named Arguments

## API


### src.embed.get_hamming_distances(genomes)
Calculate pairwise Hamming distances between the given list of genomes
and return the nonredundant array of values for use with scipy’s squareform function.
Bases other than standard nucleotides (A, T, C, G) are ignored.


* **Parameters**

    **genomes** (*list*) – a list of strings corresponding to genomes that should be compared



* **Returns**

    a list of distinct Hamming distances as a vector-form distance vector



* **Return type**

    list


```python
>>> genomes = ["ATGCT", "ATGCT", "ACGCT"]
>>> get_hamming_distances(genomes)
[0, 1, 1]
>>> genomes = ["AT-GCT", "AT--CT", "AC--CT"]
>>> get_hamming_distances(genomes)
[0, 1, 1]
```


### src.embed.make_parser()
