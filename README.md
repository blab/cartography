# Genetic cartography reveals ancestral relationships of human pathogenic viruses

**Sravani Nanduri<sup>1</sup>, Allison Black<sup>2</sup>, Trevor Bedford<sup>2,3</sup>, John Huddleston<sup>2,4</sup>**

1. Paul G. Allen School of Computer Science and Engineering, University of Washington, Seattle, WA, USA
1. Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA
1. Howard Hughes Medical Institute, Seattle, WA, USA
1. Corresponding author (jhuddles@fredhutch.org)

## Abstract

Public health studies commonly infer phylogenies from viral genome sequences to understand transmission dynamics and identify clusters of genetically-related samples.
However, viruses that reassort or recombine violate phylogenetic assumptions and require more sophisticated methods.
Even when phylogenies are appropriate, they can be unnecessary.
For example, pairwise distances between sequences can be enough to identify clusters of related samples or assign new samples to existing phylogenetic clusters.
In this work, we tested whether dimensionality reduction methods could capture known genetic groups within two human pathogenic viruses that cause substantial human morbidity and mortality and frequently reassort or recombine, respectively: seasonal influenza A/H3N2 and SARS-CoV-2.
We applied principal component analysis (PCA), multidimensional scaling (MDS), t-distributed stochastic neighbor embedding (t-SNE), and uniform manifold approximation and projection (UMAP) to sequences with well-defined phylogenetic clades and either reassortment (H3N2) or recombination (SARS-CoV-2).
For each low-dimensional embedding of sequences, we calculated the correlation between pairwise genetic and Euclidean distances in the embedding and applied a hierarchical clustering method to identify clusters in the embedding.
We measured the accuracy of these clusters compared to previously defined phylogenetic clades, reassortment clusters, or recombinant lineages.
We found that MDS maintained the strongest correlation between pairwise genetic and Euclidean distances between sequences, best captured the intermediate placement of recombinant lineages between parental lineages, and most accurately identified reassortment groups.
Clusters from t-SNE most accurately recapitulated known phylogenetic clades.
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
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov/2020-2022?m=div) colored by Nextstrain clade
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov/2020-2022?c=Nextclade_pango_collapsed&m=div) colored by collapsed Nextclade pango lineage
  - [PCA embedding](https://nextstrain.org/groups/blab/cartography/ncov/2020-2022?l=scatter&m=div&scatterX=pca1&scatterY=pca2)
  - [MDS embedding (1 and 2)](https://nextstrain.org/groups/blab/cartography/ncov/2020-2022?l=scatter&m=div&scatterX=mds1&scatterY=mds2)
  - [MDS embedding (2 and 3)](https://nextstrain.org/groups/blab/cartography/ncov/2020-2022?l=scatter&m=div&scatterX=mds2&scatterY=mds3)
  - [t-SNE embedding](https://nextstrain.org/groups/blab/cartography/ncov/2020-2022?l=scatter&m=div&scatterX=tsne_x&scatterY=tsne_y)
  - [UMAP embedding](https://nextstrain.org/groups/blab/cartography/ncov/2020-2022?l=scatter&m=div&scatterX=umap_x&scatterY=umap_y)
- Late SARS-CoV-2 (2022-2023)
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov/2022-2023?m=div) colored by Nextstrain clade
  - [Phylogeny](https://nextstrain.org/groups/blab/cartography/ncov/2022-2023?c=Nextclade_pango_collapsed&m=div) colored by collapsed Nextclade pango lineage
  - [PCA embedding](https://nextstrain.org/groups/blab/cartography/ncov/2022-2023?l=scatter&m=div&scatterX=pca1&scatterY=pca2)
  - [MDS embedding (1 and 2)](https://nextstrain.org/groups/blab/cartography/ncov/2022-2023?l=scatter&m=div&scatterX=mds1&scatterY=mds2)
  - [MDS embedding (2 and 3)](https://nextstrain.org/groups/blab/cartography/ncov/2022-2023?l=scatter&m=div&scatterX=mds2&scatterY=mds3)
  - [t-SNE embedding](https://nextstrain.org/groups/blab/cartography/ncov/2022-2023?l=scatter&m=div&scatterX=tsne_x&scatterY=tsne_y)
  - [UMAP embedding](https://nextstrain.org/groups/blab/cartography/ncov/2022-2023?l=scatter&m=div&scatterX=umap_x&scatterY=umap_y)

## Interactive figures

 - [Fig 2. **Phylogeny of early (2016–2018) influenza H3N2 HA sequences plotted by nucleotide substitutions per site on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).** Tips in the tree and embeddings are colored by their Nextstrain clade assignment.](https://blab.github.io/cartography/figures/flu-2016-2018-ha-embeddings-by-clade.html)
 - [**MDS embeddings for late (2018–2020) influenza H3N2 HA sequences showing all three components.**](https://blab.github.io/cartography/figures/flu-2016-2018-mds-by-clade.html)
 - [Fig 4. **Phylogenetic trees (left) and embeddings (right) of early (2016–2018) influenza H3N2 HA sequences colored by HDBSCAN cluster.** Normalized VI values per embedding reflect the distance between clusters and known genetic groups
(Nextstrain clades).](https://blab.github.io/cartography/figures/flu-2016-2018-ha-embeddings-by-cluster.html)
 - [S5 Fig. **Phylogeny of late (2018–2020) influenza H3N2 HA sequences plotted by nucleotide substitutions per site on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).** Tips in the tree and embeddings are colored by their Nextstrain clade assignment.](https://blab.github.io/cartography/figures/flu-2018-2020-ha-embeddings-by-clade.html)
 - [S6 Fig. **MDS embeddings for late (2018–2020) influenza H3N2 HA sequences showing all three components.**](https://blab.github.io/cartography/figures/flu-2018-2020-mds-by-clade.html)
 - [Fig 5. **Phylogenetic trees (left) and embeddings (right) of late (2018–2020) H3N2 HA sequences colored by HDBSCAN cluster.** Normalized VI values per embedding reflect the distance between clusters and known genetic groups (Nextstrain clades).](https://blab.github.io/cartography/figures/flu-2018-2020-ha-embeddings-by-cluster.html)
 - [Fig 6. **Phylogeny of early (2016–2018) influenza H3N2 HA sequences plotted by nucleotide substitutions per site on the x-axis (top) and low-dimensional embeddings of the same HA sequences concatenated with matching NA sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).** Tips in the tree and embeddings are colored by their TreeKnit Maximally Compatible Clades (MCCs) label which represents putative HA/NA reassortment groups. The first normalized VI values per embedding reflect the distance between HA/NA clusters and known genetic groups (MCCs). VI values in parentheses reflect the distance between HA-only clusters and known genetic groups.](https://blab.github.io/cartography/figures/flu-2016-2018-ha-na-embeddings-by-mcc.html)
 - [S8 Fig. **Embeddings influenza H3N2 HA-only (left) and combined HA/NA (right) showing the effects of additional NA genetic information on the
placement of reassortment events detected by TreeKnit (MCCs).**](https://blab.github.io/cartography/figures/flu-2016-2018-ha-na-all-embeddings-by-mcc.html)
 - [S9 Fig. **PCA embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/figures/flu-2016-2018-ha-na-pca-by-cluster.html)
 - [S10 Fig. **MDS embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/figures/flu-2016-2018-ha-na-mds-by-cluster.html)
 - [S11 Fig. **t-SNE embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/figures/flu-2016-2018-ha-na-tsne-by-cluster.html)
 - [S12 Fig. **UMAP embeddings for influenza H3N2 HA sequences only (top row) and HA/NA sequences combined (bottom row) showing the HA trees colored by clusters identified in each embedding (left) and the corresponding embeddings colored by cluster (right).**](https://blab.github.io/cartography/figures/flu-2016-2018-ha-na-umap-by-cluster.html)
 - [Fig 7. **Phylogeny of early (2020–2022) SARS-CoV-2 sequences plotted by number of nucleotide substitutions from the most recent common ancestor on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).** Tips in the tree and embeddings are colored by their Nextstrain clade assignment.](https://blab.github.io/cartography/figures/sarscov2-embeddings-by-Nextstrain_clade-clade.html)
 - [S13 Fig. **MDS embeddings for early SARS-CoV-2 sequences showing all three components.**](https://blab.github.io/cartography/figures/sarscov2-mds-by-Nextstrain_clade-clade.html)
 - [S15 Fig. **Phylogeny of early (2020–2022) SARS-CoV-2 sequences plotted by number of nucleotide substitutions from the most recent common ancestor on the x-axis (top) and low-dimensional embeddings of the same sequences by PCA (middle left), MDS (middle right), t-SNE (bottom left), and UMAP (bottom right).** Tips in the tree and embeddings are colored by their collapsed Nextclade pango lineage assignment.](https://blab.github.io/cartography/figures/sarscov2-embeddings-by-Nextclade_pango_collapsed-clade.html)
 - [Fig 9. **Phylogenetic trees (left) and embeddings (right) of early (2020–2022) SARS-CoV-2 sequences colored by HDBSCAN cluster. Normalized VI values per embedding reflect the distance between clusters and known genetic groups
(Nextstrain clades).**](https://blab.github.io/cartography/figures/sarscov2-embeddings-by-cluster-vs-Nextstrain_clade.html)
 - [S17 Fig. **Phylogenetic trees (left) and embeddings (right) of early (2020–2022) SARS-CoV-2 sequences colored by HDBSCAN cluster.** Normalized VI values per embedding reflect the distance between clusters and known genetic groups (collapsed Nextclade pango lineages).](https://blab.github.io/cartography/figures/sarscov2-embeddings-by-cluster-vs-Nextclade_pango_collapsed.html)
 - [Fig 10. **Phylogenetic trees (left) and embeddings (right) of late (2022–2023) SARS-CoV-2 sequences colored by HDBSCAN cluster.** Normalized VI values per embedding reflect the distance between clusters and known genetic groups
(Nextstrain clades).](https://blab.github.io/cartography/figures/sarscov2-test-embeddings-by-cluster-vs-Nextstrain_clade.html)
 - [S19 Fig. **Phylogenetic trees (left) and embeddings (right) of late (2022–2023) SARS-CoV-2 sequences colored by HDBSCAN cluster.** Normalized VI values per embedding reflect the distance between clusters and known genetic groups (collapsed Nextclade pango lineages).](https://blab.github.io/cartography/figures/sarscov2-test-embeddings-by-cluster-vs-Nextclade_pango_collapsed.html)

## Supplemental tables

 - [S1 Table. Mutations observed per embedding cluster relative to a reference genome sequence for each pathogen.](https://raw.githubusercontent.com/blab/cartography/master/docs/tables/mutation_table.csv) Each row reflects the alternate allele identified at a specific position of the given pathogen genome or gene sequence, the pathogen dataset, the embedding method, the number of clusters in the embedding with the observed mutation, and the list of distinct cluster labels with the mutation. Mutations must have occurred in at least 10 samples of the given dataset with an allele frequency of at least 50%.
 - [S2 Table. Average Euclidean distances between each known recombinant, _X_, and its parental lineages _A_ and _B_ per embedding method.](https://raw.githubusercontent.com/blab/cartography/master/sars-cov-2-nextstrain-2022-2023/results/recombinant_distances.csv) Distances include average pairwise comparisons between _A_ and _B_, _A_ and _X_, and _B_ and _X_. Additional columns indicate whether each recombinant lineage maps closer to both parental lineages (or at least one) than those parents map to each other.
 - [S3 Table. Accessions and authors from originating and submitting laboratories of seasonal influenza and SARS-CoV-2 sequences from INSDC databases.](https://raw.githubusercontent.com/blab/cartography/master/docs/tables/accessions_and_authors.tsv)

## Full analysis

## Notes for Windows  
If you are a windows user, you're going to need Linux WSL to run this. All directions on how to do this are [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

Once you have this running, you'll need to install "unzip" through linux. This is used to download and unzip the data for the MERS build. The command is below:

``` sudo apt install unzip ```

Windows Users: you CANNOT put this github repository in the Users file. Snakemake sees /U as a unicodeescape error and will not run, so please make a folder outside of the Users folder (ex. directly in the C drive) where you install this github repository, anaconda, and all other dependencies.

### Running the Cartography YAML File  
All these dependencies for this project can be installed via Conda (install conda [here](https://docs.conda.io/en/latest/miniconda.html) ) by running:

``` conda env create -f cartography.yml ```
To check the environment was created, run

``` conda-env list ```
If "cartography" is listed as one of the options, you've created the environment successfully! The asterisk will show you which environment you're currently running.

If you're running the snakemake build, there's no need to activate the environment - it does it for you. However, if you're running the evironment on its own, activate the new environment by running:

``` conda activate cartography ```

Install Javascript packages required for saving images from notebooks.

```bash conda install -c conda-forge nodejs
npm install -g vega-cli vega-lite canvas
```

### Creating a Local Repository of Cartography  
In order to work with this data and create the paper, you'll need to clone the data.
You'll first need github for your desktop. Go [Here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) for instructions on how to install github on your computer. Once github has been installed, click on the "Clone" button on the top of this repository in green, and click the clipboard icon. This allows you to copy the github repository link to your computer. You can also use the command I pasted below.

![](github_cloning_info.png)

Go back to your prompt, and navigate to a good repository to put your files (remember, "cd" allows you to type a path in, and "../" allows you to step back a directory).

Once again, do not put the repository in your "Users" folder (or any folder starting with U) - it will cause a Unicodeescape error while running the snakemake file. Finally, type

``` git clone https://github.com/blab/cartography.git ```
and wait for the repository to be created. Check to make sure everything worked by typing

``` cd cartography ```

and making sure the file location exists. Congratulations on creating a local version of the repository!

### Getting Data from NCBI for Flu  
This data is within the build itself, but if you would like to see the data we used, go to the link posted [here](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?cdate_has_day=true&cdate_has_month=true&cmd=show_query&collapse=on&country=any&fyear=2018&go=database&host=Human&lab=exclude&lineage=include&niaid=include&qcollapse=on&searchin=strain&segment=4&sequence=N&showfilters=true&sonly=on&subtype_h=3&subtype_mix=include&subtype_n=2&swine=include&tyear=2020&type=a&vac_strain=include) . This link will give you the parameters needed to get the correct data from the NCBI Influenza Database. Click the "Customize Fasta Defline" button next to the download and input


``` >{strain}|{year}-{month}-{day}|{accession}|{country}|{region} ```

It should look like this:

![](NCBI_instructions.png)

After customizing the defline, click "download", and download the data into the "data" folder of the "seasonal-flu-nextstrain" folder in your local repository of Cartography. Name the file "ncbi-h3n2-ha.fa".

Thats it! All the data has been downloaded.

### Running the Snakemake Files  
You need to install mamba now for faster downloading and updating of remote packages. Go to [this link](https://github.com/mamba-org/mamba) for full directions to do so. You will also need to install snakemake [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), which details how to install Snakemake with and without mamba (use the option for mamba)

There are two levels of snakefiles in this build. The first level lives within the respective pathogen directories (zika, mers, flu) and creates the documents, tables, and analysis needed for the paper. The second level snakefile puts the charts together with the written paper using pandoc. In order to run each part, you'll need to navigate into 1) mers-nextstrain 2) zika-nextstrain 3) seasonal-flu-nextstrain (in no particular order) and type

``` snakemake --cores 4 ```

Let the snakemake files run! They're pretty memory intensive and take a bit of time, so give the scripts time to completely run.

Once they all finish running to completion, navigate into the "docs" folder, and  run the same command (for the docs snakefile to create the paper as both a pdf and html file. All the graphs, charts, and papers made by the separate pathogen builds will also be in the docs folder. You're done!

# Running the build manually  
If you'd like to run separate parts of the build, you can definitely do so.

Creating the paper itself is done using LaTex. Simply navigate in your shell to the "docs" folder, and type
``` ./build.sh ``` for the pdf.

Running the builds for the trees can also be done separately.

For zika, navigate into the "zika-nextstrain" folder, and type
``` snakemake --use-conda --conda-frontend mamba --cores 4 ```
and let the build run. The JSON tree will be in the a top level "auspice" directory named "cartography_zika.json".

For H3N2 Ha influenza, navigate into the "seasonal-flu-nextstrain" directory in your shell, and type

``` snakemake --use-conda --conda-frontend mamba --cores 4 ```
and let the build run. The JSON tree will be in the top level "auspice" directory named "cartography_flu-seasonal-h3n2-ha-2016-2018.json".

# Build documentation  
Build the [Documentation](https://blab.github.io/cartography/):

``` bash make -C source/docs html ```

Clean the docs.

``` bash make -C source/docs clean ```

# Releasing a new version

### Information about each file

#### README.md

contains the description of the package pathogen-embed.

#### setup.py

Gives PyPi the instructions about where to find dependent packages, the authors and relevant links, etc. Also gives the entry points for the console script, which tells Pypi to call the main function of __main__.py. 

#### __init__.py

Initializes the package, creates the parser to parse the command line arguments and pass them into the embed.py function.

#### __main__.py

Calls the "run" function in __init__.py, which calls embed.py. 

#### embed.py

The main code for the package.

# To create new version 

Navigate into "source" and run 

``` python3 -m build ``` 

This creates the dist folder that gets uploaded to pypi.

``` python3 -m twine upload dist/* ```

Input the username and password, upload new dist files to pypi. Make sure the version of the dist folders does not already exist within pypi. 


