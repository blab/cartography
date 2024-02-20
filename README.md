# Dimensionality reduction distills complex evolutionary relationships in seasonal influenza and SARS-CoV-2

**Sravani Nanduri<sup>1</sup>, Allison Black<sup>2</sup>, Trevor Bedford<sup>2,3</sup>, John Huddleston<sup>2,4</sup>**

1. Paul G. Allen School of Computer Science and Engineering, University of Washington, Seattle, WA, USA
1. Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA
1. Howard Hughes Medical Institute, Seattle, WA, USA
1. Corresponding author (jhuddles@fredhutch.org)

Preprint: https://doi.org/10.1101/2024.02.07.579374

## Abstract

Public health researchers and practitioners commonly infer phylogenies from viral genome sequences to understand transmission dynamics and identify clusters of genetically-related samples.
However, viruses that reassort or recombine violate phylogenetic assumptions and require more sophisticated methods.
Even when phylogenies are appropriate, they can be unnecessary or difficult to interpret without specialty knowledge.
For example, pairwise distances between sequences can be enough to identify clusters of related samples or assign new samples to existing phylogenetic clusters.
In this work, we tested whether dimensionality reduction methods could capture known genetic groups within two human pathogenic viruses that cause substantial human morbidity and mortality and frequently reassort or recombine, respectively: seasonal influenza A/H3N2 and SARS-CoV-2.
We applied principal component analysis (PCA), multidimensional scaling (MDS), t-distributed stochastic neighbor embedding (t-SNE), and uniform manifold approximation and projection (UMAP) to sequences with well-defined phylogenetic clades and either reassortment (H3N2) or recombination (SARS-CoV-2).
For each low-dimensional embedding of sequences, we calculated the correlation between pairwise genetic and Euclidean distances in the embedding and applied a hierarchical clustering method to identify clusters in the embedding.
We measured the accuracy of clusters compared to previously defined phylogenetic clades, reassortment clusters, or recombinant lineages.
We found that MDS maintained the strongest correlation between pairwise genetic and Euclidean distances between sequences and best captured the intermediate placement of recombinant lineages between parental lineages
Clusters from t-SNE most accurately recapitulated known phylogenetic clades and recombinant lineages.
Both MDS and t-SNE accurately identified reassortment groups.
We show that simple statistical methods without a biological model can accurately represent known genetic relationships for relevant human pathogenic viruses.
Our open source implementation of these methods for analysis of viral genome sequences can be easily applied when phylogenetic methods are either unnecessary or inappropriate.

## Phylogenetic trees and embeddings

Explore the phylogenetic trees and embeddings on Nextstrain.

- Early influenza H3N2 HA (2016-2018)
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018?m=div) colored by Nextstrain clade
  - [PCA embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018?l=scatter&m=div&scatterX=pca1&scatterY=pca2)
  - [MDS embedding (1 and 2)](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018?l=scatter&m=div&scatterX=mds1&scatterY=mds2)
  - [MDS embedding (2 and 3)](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018?l=scatter&m=div&scatterX=mds2&scatterY=mds3)
  - [t-SNE embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018?l=scatter&m=div&scatterX=tsne_x&scatterY=tsne_y)
  - [UMAP embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018?l=scatter&m=div&scatterX=umap_x&scatterY=umap_y)
- Late influenza H3N2 HA (2018-2020)
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2018-2020?m=div) colored by Nextstrain clade
  - [PCA embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2018-2020?l=scatter&m=div&scatterX=pca1&scatterY=pca2)
  - [MDS embedding (1 and 2)](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2018-2020?l=scatter&m=div&scatterX=mds1&scatterY=mds2)
  - [MDS embedding (2 and 3)](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2018-2020?l=scatter&m=div&scatterX=mds2&scatterY=mds3)
  - [t-SNE embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2018-2020?l=scatter&m=div&scatterX=tsne_x&scatterY=tsne_y)
  - [UMAP embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2018-2020?l=scatter&m=div&scatterX=umap_x&scatterY=umap_y)
- Influenza H3N2 HA and NA (2016-2018)
  - [HA phylogeny](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018-reassort) colored by maximum compatibility clades (MCCs) representing HA/NA reassortment groups
  - [NA phylogeny](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-na-2016-2018-reassort) colored by MCCs
  - [HA/NA tangletree](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018-reassort:groups/blab/cartography/flu-seasonal-h3n2-na-2016-2018-reassort) colored by MCCs
  - [PCA embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018-reassort?l=scatter&scatterX=pca1&scatterY=pca2)
  - [MDS embedding (1 and 2)](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018-reassort?l=scatter&scatterX=mds1&scatterY=mds2)
  - [MDS embedding (2 and 3)](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018-reassort?l=scatter&scatterX=mds2&scatterY=mds3)
  - [t-SNE embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018-reassort?l=scatter&scatterX=tsne_x&scatterY=tsne_y)
  - [UMAP embedding](https://nextstrain.org/groups/blab/cartography/flu-seasonal-h3n2-ha-2016-2018-reassort?l=scatter&scatterX=umap_x&scatterY=umap_y)
- Early SARS-CoV-2 (2020-2022)
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov-2020-2022?m=div) colored by Nextstrain clade
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov-2020-2022?c=Nextclade_pango_collapsed&m=div) colored by collapsed Nextclade pango lineage
  - [PCA embedding](https://nextstrain.org/groups/blab/cartography/ncov-2020-2022?l=scatter&m=div&scatterX=pca1&scatterY=pca2)
  - [MDS embedding (1 and 2)](https://nextstrain.org/groups/blab/cartography/ncov-2020-2022?l=scatter&m=div&scatterX=mds1&scatterY=mds2)
  - [MDS embedding (2 and 3)](https://nextstrain.org/groups/blab/cartography/ncov-2020-2022?l=scatter&m=div&scatterX=mds2&scatterY=mds3)
  - [t-SNE embedding](https://nextstrain.org/groups/blab/cartography/ncov-2020-2022?l=scatter&m=div&scatterX=tsne_x&scatterY=tsne_y)
  - [UMAP embedding](https://nextstrain.org/groups/blab/cartography/ncov-2020-2022?l=scatter&m=div&scatterX=umap_x&scatterY=umap_y)
- Late SARS-CoV-2 (2022-2023)
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov-2022-2023?m=div) colored by Nextstrain clade
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov-2022-2023?c=Nextclade_pango_collapsed&m=div) colored by collapsed Nextclade pango lineage
  - [PCA embedding](https://nextstrain.org/groups/blab/cartography/ncov-2022-2023?l=scatter&m=div&scatterX=pca1&scatterY=pca2)
  - [MDS embedding (1 and 2)](https://nextstrain.org/groups/blab/cartography/ncov-2022-2023?l=scatter&m=div&scatterX=mds1&scatterY=mds2)
  - [MDS embedding (2 and 3)](https://nextstrain.org/groups/blab/cartography/ncov-2022-2023?l=scatter&m=div&scatterX=mds2&scatterY=mds3)
  - [t-SNE embedding](https://nextstrain.org/groups/blab/cartography/ncov-2022-2023?l=scatter&m=div&scatterX=tsne_x&scatterY=tsne_y)
  - [UMAP embedding](https://nextstrain.org/groups/blab/cartography/ncov-2022-2023?l=scatter&m=div&scatterX=umap_x&scatterY=umap_y)

## Interactive figures

### Main figures

 - [Fig 2. **Phylogeny of early (2016–2018) influenza H3N2 HA sequences plotted by nucleotide substitutions per site on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).**](https://blab.github.io/cartography/flu-2016-2018-ha-embeddings-by-clade.html) Tips in the tree and embeddings are colored by their Nextstrain clade assignment. Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line colors represent the clade membership of the most ancestral node in the pair of nodes connected by the segment. Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [Fig 4. **Phylogenetic trees (left) and embeddings (right) of early (2016–2018) influenza H3N2 HA sequences colored by HDBSCAN cluster.**](https://blab.github.io/cartography/flu-2016-2018-ha-embeddings-by-cluster.html) Normalized VI values per embedding reflect the distance between clusters and known genetic groups (Nextstrain clades). Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [Fig 5. **Phylogenetic trees (left) and embeddings (right) of late (2018–2020) H3N2 HA sequences colored by HDBSCAN cluster.**](https://blab.github.io/cartography/flu-2018-2020-ha-embeddings-by-cluster.html) Normalized VI values per embedding reflect the distance between clusters and known genetic groups (Nextstrain clades). Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [Fig 6. **Phylogeny of early (2016–2018) influenza H3N2 HA sequences plotted by nucleotide substitutions per site on the x-axis (top) and low-dimensional embeddings of the same HA sequences concatenated with matching NA sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).**](https://blab.github.io/cartography/flu-2016-2018-ha-na-embeddings-by-mcc.html) Tips in the tree and embeddings are colored by their TreeKnit Maximally Compatible Clades (MCCs) label which represents putative HA/NA reassortment groups. The first normalized VI values per embedding reflect the distance between HA/NA clusters and known genetic groups (MCCs). VI values in parentheses reflect the distance between HA-only clusters and known genetic groups. "A2" and "A2/re" labels indicate a known reassortment event ([Potter et al. 2019](https://doi.org/10.1093/ve/vez046)).
 - [Fig 7. **Phylogeny of early (2020–2022) SARS-CoV-2 sequences plotted by number of nucleotide substitutions from the most recent common ancestor on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).**](https://blab.github.io/cartography/sarscov2-embeddings-by-Nextstrain_clade-clade.html) Tips in the tree and embeddings are colored by their Nextstrain clade assignment. Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [Fig 9. **Phylogenetic trees (left) and embeddings (right) of early (2020–2022) SARS-CoV-2 sequences colored by HDBSCAN cluster.**](https://blab.github.io/cartography/sarscov2-embeddings-by-cluster-vs-Nextstrain_clade.html) Normalized VI values per embedding reflect the distance between clusters and known genetic groups (Nextstrain clades). Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [Fig 10. **Phylogenetic trees (left) and embeddings (right) of late (2022–2023) SARS-CoV-2 sequences colored by HDBSCAN cluster.**](https://blab.github.io/cartography/sarscov2-test-embeddings-by-cluster-vs-Nextstrain_clade.html) Normalized VI values per embedding reflect the distance between clusters and known genetic groups (Nextstrain clades).

### Supplemental figures

 - [S4 Fig. **MDS embeddings for early (2016–2018) influenza H3N2 HA sequences showing all three components.**](https://blab.github.io/cartography/flu-2016-2018-mds-by-clade.html) Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line colors represent the clade membership of the most ancestral node in the pair of nodes connected by the segment. Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [S6 Fig. **Phylogeny of late (2018–2020) influenza H3N2 HA sequences plotted by nucleotide substitutions per site on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).**](https://blab.github.io/cartography/flu-2018-2020-ha-embeddings-by-clade.html) Tips in the tree and embeddings are colored by their Nextstrain clade assignment. Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line colors represent the clade membership of the most ancestral node in the pair of nodes connected by the segment. Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [S7 Fig. **MDS embeddings for late (2018–2020) influenza H3N2 HA sequences showing all three components.**](https://blab.github.io/cartography/flu-2018-2020-mds-by-clade.html) Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line colors represent the clade membership of the most ancestral node in the pair of nodes connected by the segment. Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [S9 Fig. **Embeddings influenza H3N2 HA-only (left) and combined HA/NA (right) showing the effects of additional NA genetic information on the
placement of reassortment events detected by TreeKnit (MCCs).**](https://blab.github.io/cartography/flu-2016-2018-ha-na-all-embeddings-by-mcc.html)
 - [S10 Fig. **PCA embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/flu-2016-2018-ha-na-pca-by-cluster.html)
 - [S11 Fig. **MDS embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/flu-2016-2018-ha-na-mds-by-cluster.html)
 - [S12 Fig. **t-SNE embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/flu-2016-2018-ha-na-tsne-by-cluster.html)
 - [S13 Fig. **UMAP embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/flu-2016-2018-ha-na-umap-by-cluster.html)
 - [S14 Fig. **MDS embeddings for early SARS-CoV-2 sequences showing all three components.**](https://blab.github.io/cartography/sarscov2-mds-by-Nextstrain_clade-clade.html) Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [S15 Fig. **Phylogeny of early (2020–2022) SARS-CoV-2 sequences plotted by number of nucleotide substitutions from the most recent common ancestor on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).**](https://blab.github.io/cartography/sarscov2-embeddings-by-Nextclade_pango_collapsed-clade.html) Tips in the tree and embeddings are colored by their collapsed Nextclade pango lineage assignment. Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [S17 Fig. **Phylogenetic trees (left) and embeddings (right) of early (2020–2022) SARS-CoV-2 sequences colored by HDBSCAN cluster.**](https://blab.github.io/cartography/sarscov2-embeddings-by-cluster-vs-Nextclade_pango_collapsed.html) Normalized VI values per embedding reflect the distance between clusters and known genetic groups (collapsed Nextclade pango lineages). Line segments in each embedding reflect phylogenetic relationships with internal node positions calculated from the mean positions of their immediate descendants in each dimension (see Methods). Line thickness scales by the square root of the number of leaves descending from a given node in the phylogeny.
 - [S19 Fig. **Phylogenetic trees (left) and embeddings (right) of late (2022–2023) SARS-CoV-2 sequences colored by HDBSCAN cluster.**](https://blab.github.io/cartography/sarscov2-test-embeddings-by-cluster-vs-Nextclade_pango_collapsed.html) Normalized VI values per embedding reflect the distance between clusters and known genetic groups (collapsed Nextclade pango lineages).

## Supplemental tables

 - [S1 Table. Optimal cluster thresholds per pathogen, known genetic group type, and embedding method based on normalized variation of information
(VI) distances calculated from early pathogen datasets.](https://raw.githubusercontent.com/blab/cartography/master/manuscript/tables/total_HDBSCAN_table.csv) Smaller VI values indicate fewer differences between HDBSCAN clusters and known genetic groups. VI of 0 indicates identical clusters and 1 indicates maximally different clusters. Threshold refers to the minimum Euclidean distance between initial clusters for HDBSCAN to consider them as distinct clusters. We apply these optimal thresholds per pathogen, known genetic group type, and method to find clusters in corresponding late datasets for each pathogen.
 - [S2 Table. Number of clusters, transitions between clusters in the phylogeny, and excess transitions indicating non-monophyletic groups per pathogen and embedding.](https://raw.githubusercontent.com/blab/cartography/master/manuscript/tables/monophyletic_clusters.csv) Embeddings without any excess transitions reflect monophyletic groups in the corresponding pathogen phylogeny.
 - [S3 Table. Mutations observed per embedding cluster relative to a reference genome sequence for each pathogen.](https://raw.githubusercontent.com/blab/cartography/master/manuscript/tables/mutation_table.csv) Each row reflects the alternate allele identified at a specific position of the given pathogen genome or gene sequence, the pathogen dataset, the embedding method, the number of clusters in the embedding with the observed mutation, and the list of distinct cluster labels with the mutation. Mutations must have occurred in at least 10 samples of the given dataset with an allele frequency of at least 50%.
 - [S4 Table. Average Euclidean distances between each known recombinant, _X_, and its parental lineages _A_ and _B_ per embedding method.](https://raw.githubusercontent.com/blab/cartography/master/manuscript/tables/sars-cov-2-recombinant-distances.csv) Distances include average pairwise comparisons between _A_ and _B_, _A_ and _X_, and _B_ and _X_. Additional columns indicate whether each recombinant lineage maps closer to both parental lineages (or at least one) than those parents map to each other.
 - [S5 Table. Accessions and authors from originating and submitting laboratories of seasonal influenza and SARS-CoV-2 sequences from INSDC databases.](https://raw.githubusercontent.com/blab/cartography/master/manuscript/tables/accessions_and_authors.tsv)

## Full analysis

### Installation

First, [install Conda with the Miniconda distribution](https://docs.conda.io/en/latest/miniconda.html).
Until Bioconda supports modern Mac CPUs, Mac users with M1/M2 CPUs (the ARM64 architecture) need to install the Mac Intel x86 Miniconda distribution and [install Rosetta](https://support.apple.com/en-us/HT211861), so the workflow can run under Mac's emulation mode.

After installing Conda, create the environment for this project.

```bash
conda env create -f cartography.yml
```

Activate the environment prior to running the workflow below.

```bash
conda activate cartography
```

Next, you need to [install Julia](https://julialang.org/downloads/) and then [install TreeKnit](https://pierrebarrat.github.io/TreeKnit.jl/) following the instructions to install the "CLI" version.
The TreeKnit binary installs in your home directory, by default, in the path `~/.julia/bin/treeknit`.
This path is what the project's workflow calls to run TreeKnit.

### Notes for Windows users

If you are a Windows user, you will need to [install WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10) to run this project's workflow.
You _cannot_ put this github repository in the Users file. Snakemake sees /U as a unicodeescape error and will not run, so please make a folder outside of the Users folder (ex. directly in the C drive) where you install this github repository, anaconda, and all other dependencies.

### Run the full analysis

Run the full analysis for the project which includes simulations, analysis of natural populations, and generation of the manuscript and its figures and tables.

```
snakemake \
    --use-conda \
    --conda-frontend conda \
    --cores all
```

This is a complex workflow, so it will take several hours to run.
