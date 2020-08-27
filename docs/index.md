---
title: "Cartography"
author: "Sravani Nanduri"
---


# Links to Graphs/Charts Cartography

## Zika Charts
[MDS Clickable](https://blab.github.io/cartography/MDS.html)\
[PCA Brush](https://blab.github.io/cartography/PCABrush.html)\
[TSNE Clickable](https://blab.github.io/cartography/TSNEClickable.html)\
[TSNE Learning Rates](https://blab.github.io/cartography/TSNELearningRates.html)\
[UMAP Clickable](https://blab.github.io/cartography/UMAPClickable.html)\
[PCA Scaled Clickable](https://blab.github.io/cartography/PCAScaledClickable.html)\

## Flu Charts
[MDS Brush](https://blab.github.io/cartography/MDSFluBrush.html)\
[TSNE Clickable](https://blab.github.io/cartography/TSNEFluClickable.html)\
[PCA Brush](https://blab.github.io/cartography/PCAFluBrush.html)\
[PCA Clickable](https://blab.github.io/cartography/PCAFluClickable.html)\
[UMAP Brush](https://blab.github.io/cartography/UMAPFluBrush.html)\

# DESCRIPTION OF THE PROBLEM:

Phylogenetic inference is a fundamental tool for understanding genealogical relationships among human pathogenic viruses. 
However, recombination and reassortment of viral genomes invalidates basic phylogenetic assumptions of inheritance and requires more sophisticated approaches.
One approach is to split a genome into multiple phylogenies to model the evolution of the nonrecombinant fragments.
This is done using a genetic algorithm that scans strains for recombination breakpoints, quantifies and analyzes the impact of recombination at each one, and splits the phylogeny at its most important breakpoints [@kosakovsky_pond_posada_gravenor_woelk_frost_2006].
Finding recombination breakpoints relies on the detection of a recombination signal through methods such as CHIMAERA and LARD.
Both CHIMAERA and LARD use split decomposition, a method which depicts parallel edges between sequences if there are conflicting phylogenetic signals in the data [@Posada13757][@martin_murrell_khoosal_muhire_2017].
An alternate strategy is to compare viral genomes with alternative methods that do not make the same strong assumptions as phylogenetic inference (e.g., PCA, MDS, etc.).
PCA has been used to estimate and model ancestry in unrelated individuals and plot peoples genomes to reveal patterns in national origin [@alexander_novembre_lange_2009][@GenesMirrorGeography2008], and was also used to genotype major classes of structural variants - structural differences such as deletion, duplication, and novel sequence insertion - in diverse populations to map their population stratification [@sudmant_1000_gene_paper].
PCA was also used to reveal Zika’s genetic diversity and spread in the Americas by assessing clustering of multidimensional genetic data [@H.C.109348].
MDS has been applied to h3n2 sequences to inspect relationships between all gene segments, which is closely related to the subject of this paper, for the difference that Rambaut et .al. 2008 looks at between-gene diversity rather than within-gene [@rambaut_pybus_nelson_viboud_taubenberger_holmes_2008].
PCA, t-SNE, and UMAP have all been used to better capture both discrete and continuous patterns of variation in human genomes across a genetic continuum, and the embeddings were able to show relationships between genotype, phenotype, and geography [diaz-papkovich_anderson-trocmé_ben-eghan_gravel_2019]. 
While Diaz-Papkovich et.al. and Metsky et.al. explored qualitative measurements of embedding accuracy and fitness, this paper will go beyond that by establishing quantitative measurements as to the fit and accuracy of the embeddings to further bridge the gap between visualization and statistical testing.  
This paper will also give insight into different reduction techniques, and will discuss both their limitations and strengths in the realm of viral data. 
We present a novel approach to understanding relationships among viral genomes by transforming genomic data and then using dimensionality reduction methods such as PCA, MDS, t-SNE, and UMAP. 
We investigate the degree to which this method can recapitulate known phylogenetic relationships for viruses whose genomes are phylogenetically tractable (we used influenza h3n2 and Zika).
We apply this method to viruses whose genomes are known to undergo substantial recombination, such as MERS and SARS-CoV-2 to assess how well each method is able to reconstruct previously identified biologically-meaningful clusters.

Recombination: occurs when at least two viral genomes co-infect the same host cell and exchange genetic segments. 
Shuffling/reassortment, a particular type of recombination, occurs in viruses with segmented genomes, which by interchanging complete genome segments, gives rise to new segment combinations [@pérez-losada_arenas_galán_palero_gonzález-candelas_2014].


# METHODS:

### Materials:

The analysis environment can be recreated using conda and all installation instructions are available on [Cartography’s github](https://github.com/blab/cartography) .

### Methods:

The genome data we used for h3n2 HA influenza is from the NCBI Influenza database. 
We used [this search](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?cdate_has_day=true&cdate_has_month=true&cmd=show_query&collapse=on&country=any&defline_saved=%3E%7Baccession%7D%20%7Bstrain%7D%20%7Byear%7D/%7Bmonth%7D/%7Bday%7D%20%7Bsegname%7D&fyear=2015&go=database&host=Human&lab=exclude&lineage=include&niaid=include&qcollapse=on&searchin=strain&segment=4&sequence=N&showfilters=true&sonly=on&subtype_h=3&subtype_mix=include&subtype_n=2&swine=include&tyear=2020&type=a&vac_strain=include). Clades were defined by reasonable phylogenetic signal. 
The Zika data was curated by Allison Black, with sequences from Genbank and the Bedford Lab. Clades were defined by regionally important introductions as well as by reasonable phylogenetic signal in terms of mutations on branches.  
We analyzed Influenza A/h3n2, Zika, and MERS genomes and created a FASTA file of multiple sequence alignments with MAFFT v7.407 [@Katoh2002] via augur align [@Hadfield2018] and phylogenies with IQ-TREE v1.6.10 [@Nguyen2014] via augur tree version 5.4.1.

We used two different methods of transforming the data; Scaling and centering the data, and a Hamming distance similarity matrix.
For Scaling and Centering the data, we performed PCA on the matrix of nucleotides from the multiple sequence alignment using scikit-learn [@jolliffe_cadima_2016]. 
An explained variance plot was created to determine the amount of PCs created, which is in the supplementary figures section.
We also ran UMAP on the first 10 leading PCs to see if patterns not fully captured by the embedding were seen through UMAP clustering. 
The plots and analysis for UMAP on the first 10 PCs is available in the supplemental figures section.

For Hamming distance, we created a similarity matrix. 
By comparing every genome with every other genome and clustering based on their Hamming distance, distance-based methods take the overall structure of the multidimensional data and groups together genomes that have similar differences.
This means the data is clustered by genetic diversity (in a phylogenetic tree genetic diversity is categorized using clades).
Each genome was split into separate nucleotides and compared with other nucleotides in the same site on other genomes.
We only counted a difference between the main nucleotide pairs (AGCT) -- gaps (N) were not.
This is because some sequences were significantly shorter than others, and a shorter strain does not necessarily mean complete genetic dissimilarity, which is what counting gaps implied. 

We reduced the similarity distance matrix through MDS, t-SNE, and UMAP, plotted using [Altair](https://altair-viz.github.io/) ,and colored by clade assignment.
Clade membership metadata was provided by a .json build of the influenza h3n2 tree (the build can be found at https://github.com/blab/cartography/tree/master/notebooks/Data).
The 3 different dimensionality reduction techniques are ordered below by publication date: 
- [MDS](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.MDS.html)
- [t-SNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)
- [UMAP](https://umap-learn.readthedocs.io/en/latest/)

The plots of the full 10 PCs for PCA and the first 6 components for MDS are available in the supplemental figures section.

We tested different learning rates and perplexity values for t-SNE and UMAP by plugging in a variety of permutations of parameters to determine the best embedding. 
These plots are in the supplemental figures section. 
The default value of a 200.0 learning rate and a 30 perplexity provided good results, but we ended up going with a 25.95 perplexity due to various computational advantages. 
Similarly, we tested different parameter choices for UMAP, with the clearest results generated by specifying 200 nearest neighbours and a “minimum distance” between points in low dimensions of 0.05. 
While tuning these parameters will not change qualitative results, it can help make patterns easier to identify. For example, the more nearest neighbors, the higher the computational load, and while smaller minimum distances can break connectivity between clusters, they will not change the groupings of individuals.

# Tuning of UMAP: goes here once done


To further analyze these embeddings’ ability to accurately capture the multidimensional data, we made two separate plots: hamming vs euclidean distance scatterplots with a LOESS best fit line, and within vs between clade violin plots per embedding.

Hamming distance vs euclidean distance scatterplots:
 
Hamming distance vs Euclidean distance plots assess the local and global structure of the embedding as well as assess the overall strength of the embedding's recapitulation.
The Hamming distance between nucleotide sequences is plotted on the x axis, and the euclidean distance between the points in the embedding are plotted on the y axis.
By plotting these distance measurements, we can observe how correlated the dataset is.
The higher the correlation, the better a function can describe the relationship between the Hamming distance value and the euclidean distance value. 
In this way, constant correlation in a plot reveals that the embedding tends to capture and retain local patterns rather than global, and a splayed structure points to global structure preservation over local. 
Therefore, the closer the Pearson Coefficient is to 1, the better the embedding is at preserving genetic dissimilarity in euclidean space.
The LOESS line drawn through the plot assesses the best fit function for the embedding.

To quantify the patters seen in the scatterplot further, we bootstrapped our scatterplot _________________________________


Between vs Within clade Violin plots:

The Between vs Within clade Violin plots visually represent how well Euclidean distances can distinguish virus genomes from different clades.
In other words, it describes the probability that a certain Euclidean distance can be used to classify a given pair of genomes as within vs between clades. 
The density of the violin plot at a specific distance gives the relative probability of two strains with that distance being in same or different clades.
The median to median ratio is an indicator of how well the embedding clusters the data.
The larger the ratio between the medians of between vs within clade violin plots, the better the embedding is at clustering and compartmentalizing data into their clades.
To create this plot, the matrix of euclidean distances for each embedding was flattened, and each comparison was labeled as a “within clade” or “between clade” comparison using the clade assignments from the .json build of the tree.
Violin plots were made using [seaborn](https://seaborn.pydata.org/) , separated by clade status and euclidean distance on the y axis. 

 To quantify the patterns seen in the plot beyond a median ratio, we ran a classification test to test the significance of differences in median for between and within clusters. ________________
 This test would answer if viral genomes from the same clades have smaller Euclidean distances in a given embedding than viral genomes from different clades. 
 It also answers if any of these embeddings classify pairs of viral genomes better than genetic Hamming distance, which helps us decide if creating the computationally intensive embeddings is worthwhile. 
 
 We constructed a classification test by ________________________________
 
 

# RESULTS:

### EXPECTATIONS FOR PCA, MDS, t-SNE, and UMAP

Principal Component Analysis (PCA) reduces multidimensional data, increasing interpretability while minimizing information loss[@jolliffe_cadima_2016] . 
PCA relies on linear assumptions, does not affect the scale of the data, and does not normalize the data as part of the algorithm.
PCA preserves long range distances but hides finer-scale details. Because PCA is almost entirely focused on retaining the global structure and variance of the data, and one of its limitations is revealing patterns locally. 
PCA is not an algorithm to be used on a similarity matrix, and is instead intended for transformed and normalized multidimensional data.
In the context of this paper, PCA will be used on transformed and normalized genetic data and not on the similarity matrix described above.

Multidimensional Scaling (MDS) refers to statistical techniques that reduce the complexity of a data set by quantifying similarity judgments, which increases the interpretability of local relational structures mired in the dataset [@hout_papesh_goldinger_2012] .
A limitation to MDS is that only one symmetric matrix is allowed as input, and the scale of measurement is non-numerical.
MDS preserves global patterns over local, but the algorithm’s importance on translating dissimilarity to distance does preserve some larger patterns in local structure as well.
In the context of this paper, MDS will cluster data into sparse “sections” of the map while not creating actual clusters. 

t-distributed Stochastic Neighbor Embedding (t-SNE) visualizes high-dimensional data by giving each datapoint a location in a 2 to 3 dimensional map.
t-SNE is focused largely on local structure over global structure, and t-SNE’s projection of clusters and distances between clusters are not analogous to dissimilarity - in other words, t-SNE focuses heavily on projecting similarity rather than dissimilarity [@maaten2008visualizing].
Because t-SNE reduces data’s dimensionality based on local properties of data, data with intrinsically high dimensional structure will not be projected accurately.
In the context of this paper, t-SNE will create tight clusters that clearly indicate genetic similarity, but will not create an accurate global picture of the data.

Uniform Manifold Approximation and Projection (UMAP) is a manifold learning technique for dimension reduction based in Riemannian geometry and algebraic topology [@lel2018umap]. 
The end result is a patchwork of low-dimensional representations of neighbourhoods that groups genetically similar strains together on a local scale while better preserving long-range topological connections to more distantly related strains.
Some limitations include its lack of maturity - this novel technique does not have firmly established or robust practices and libraries to use UMAP best.
In the context of this paper, UMAP will reveal a tightly clustered set of data that retains both the global structure of the data and the clusters and similarities present at the local level. 

# EXPECTATIONS FOR VIRUSES:

## Influenza:

h3n2 Influenza in this project is used as a proof of concept as h3n2 HA influenza only reassorts and does not recombine. 
We use h3n2's HA sequences as they have a relatively high mutation rate compared to the other gene segments, it encodes a protein that is a target of human immunity, and has traditionally been used for analysis of influenza evolution. 
As these sequences are biologically relevant, short, and do not recombine, the genomes can be reasonably assigned to phylogenetic clades. 
Therefore, it can be assumed that h3n2 HA influenza is a good test case for Cartography.

## SUMMARY OF RESULTS FOR INFLUENZA 

PCA: There were visually identifiable clusters within the data (Figure 1B), with clades A1b/135K and A3 being overlapped by “other” clades. 
Overall, the graph segmented pretty well with some clades being broken into more than one cluster, such as A1b/131K and A1b/135K.
In the euclidean and Hamming distance scatter plot for PCA (Figure 2A), the points are at a constant distance from the LOESS line as genetic distance increases. 
This correlation’s fit can be expressed through the scatterplot's Pearson Coefficient of .739. 
The LOESS line was fairly linear, which, given that PCA retains global patterns and therefore translates genetic distance as accurately to euclidean distance as possible, upholds preexisting beliefs about the algorithm.
Euclidean distance can be used with some confidence to distinguish strains from similar and different clades (Figure 3A). 
In Figure 3A, the within density in the violin plot is in lower distances than the density of the between clade violin plot, their medians have a large difference between them, and the plot resembles the same shape as the genetic violin plot. 
These observations reveal that the embedding did not change the structure of the data, but instead exaggerated distances on a global scale to make the data more visually useful. 
However, the data did not seem to show any new patterns of information and clustered the data points almost identically to the .json rendering of the tree; similar distance and placement of clade clusters corroborates this hypothesis. _______________



MDS: MDS, while there are visually identifiable clusters in the embedding, has most of its points and clusters overlapping each other in a single mass in the middle of the graph (Figure 1C). 
In particular, the clusters of clades A2/re, A3, A2, A1b/135K, and A1b/131K all overlapped, as they are the most genetically similar clades according to the tree. 
In the Euclidean and Hamming distance scatterplot (Figure 2B), the points stay relatively close to the LOESS curve, and the points correlation according to the LOESS line is corroborated by the Pearson Coefficient of .616 for the plot. 
This linear pattern in the data points is typical of an algorithm that preserves mostly global relationships while still segmenting and compartmentalizing the data better than the distance matrix does on its own. 
Euclidean distance does not do an accurate job distinguishing viruses from similar and different clades (Figure 3B). 
In Figure 3B, the within density in the violin plot is concentrated in lower distances, and this is also seen in the density of the between clade violin plot. 
Because the MDS embedding put all the data points in a central mass in the middle of the embedding, this finding is not surprising. 
Quantitatively, this inaccuracy can be described by the P-value of ___________



t-SNE: The t-SNE embedding was very accurate at clustering the data by genetic diversity. 
t-SNE also found a clade breakpoint previously unseen; one of the clades (A1b/131K) was split into three clusters that, when looking at NextStrain’s updated Ha phylogeny, slices the clade in three clean pieces. 
Looking at newly defined branches in Nextstrain, the “slicing” in the embedding revealed patterns in the tree that were not seen, and this foreshadows Cartography’s possible competency with pinpointing future clade break points.
In the euclidean and Hamming distance scatterplot for t-SNE (Figure 2C), the points splay out from the LOESS line as genetic distance increases. 
The data’s correlation to the LOESS line can be expressed through the scatterplot's Pearson Coefficient of .427. 
This increasing divergence between genetic and euclidean distance is characteristic of an embedding that exposes local patterns in the data over global structure.
Euclidean distance is a fairly strong measure for distinguishing strains as between or within clades. 
The density of the within violin plot was mostly concentrated at lower euclidean distances while the between clade violin plot’s density was at higher euclidean distances (Figure 3C). 
The medians given by the within vs between plots were also fairly distant from each other. 
The P-value corroborates this accuracy, ____________



UMAP: The UMAP embedding did a very accurate job clustering the data by genetic diversity. 
There was some overlap between A2/re and A2, and a few data points in the wrong cluster, but overall the UMAP projection densely packed together strains within a clade and put lots of distance for strains between clades. 
The data points did splay out from the LOESS line in the beginning, but far less so compared to t-SNE - the more distinguishing characteristics were the clusters shown in the scatter plot itself, revealing UMAP’s approach to visualizing local and global patterns equally (Figure 2D). 
The correlation of the data points to the LOESS line is described by UMAP’s Pearson Coefficient of .369. 
The LOESS line plateaus at a genetic distance of 30, revealing that UMAP sets a threshold for genetic diversity and keeps all points above that threshold the same distance apart in euclidean space.
Euclidean distance is a strong measure for distinguishing strains as between or within clades. 
The density of the within violin plot was almost completely concentrated at euclidean distances lower than 10 while the between clade violin plot’s density was above that distance (Figure 3D). This accuracy can we quantified through the P-value of ________________ 



## Figure One

:::{.out .html}
<iframe src="https://blab.github.io/cartography/FullLinkedChartClickableFlu.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](FullLinkedChartClickableFlu.pdf)
:::


## Figure Two

![](FullScatterplotFlu.png)
 
## Figure Three

![](FullViolinPlotFlu.png)




## Zika:

Zika: Zika in this project is used as a test case. 
While h3n2 Influenza is a globally distributed virus that has caused infections seasonally for decades, Zika is a fairly new human pathogenic virus that has a restricted geographic distribution that recapitulates the patterns of viral transmission. 
Therefore, Zika is better compartamentalized by region than by clade. 
It can be reasonably assumed that zika is a good test case for Cartography.

## SUMMARY OF RESULTS FOR ZIKA: 

PCA: The PCA plot for Zika colored by clade (Figure 4B) did not reveal any visually identifiable clusters in the data, and the data’s overall structure also did not reveal anything interesting. 
Visually, there was no relationship between Hamming Distance and Euclidean Distance within the PCA embedding, shown by the completely random placement of points on the scatterplot (Figure 5A). 
To quantify this visual observation, the Pearson Coefficient was .005, which reveals close to no correlation. 
The lack of visible clustering in this embedding reveals that scaling and centering nucleotide data does not capture the genetic diversity present between each genome.  
Euclidean distance does not help distinguish viral genomes by genetic diversity (Figure 6A). ___________________________-

MDS: The embedding seemed to differentiate between the clusters in the data on a very global scale (Figure 4C). 
There were a few visually identifiable clusters in the data, such as clusters containing clades 1, 10, and 2, but there was a single mass of points on the graph in the center where clades 3,4,5, and 6 were slightly overlapping.
There was segmentation in the colors of the data and their overlap in this single mass of points, however - clade 4 is layered over clade 5. Because MDS tends to reveal global structure over local, MDS tries to as accurately as possible preserve the data’s overall structure. 
In the Euclidean and Hamming distance scatterplot, the points begin relatively close to the LOESS curve, but begin to move farther and farther away as the genetic distance increases, and this divergence is corroborated by the Pearson Coefficient of .376 for the plot. 
This “splayed” pattern as Euclidean and Hamming distance begin to diverge increasingly is typical of an algorithm that preserves both local and global relationships - global due to the linear relationship between euclidean and hamming distance, and local from the splay pattern - which reveals that MDS seems to be segmenting and compartmentalizing the data better than the distance matrix does on its own. 
Euclidean distance does a fairly accurate job distinguishing viruses from similar and different clades (Figure 6B). 
In Figure 6B, the within density in the violin plot is in lower distances than the density of the between clade violin plot, which means their medians have a large difference between them. 
In the permutation test, we obtained a p-value of  ____________________________-


t-SNE: This embedding does a very good job of clustering the data; every clade is a visually identifiable cluster, and clades more similarly related genetically are closer together in the plot (Figure 4D). 
The data points from different clades did not overlap each other much, and the embedding did revealed different relationships than the rendering of the tree, giving less euclidean space between points that were more genetically diverse and creating a “threshold” for how much genetic diversity would impact euclidean distance. 
In the euclidean and Hamming distance scatterplot for t-SNE (Figure 5C), the points stay at a constant distance from the LOESS line as genetic distance increases. 
This correlation can be expressed through the scatterplot's Pearson Coefficient of .455. 
The LOESS line plateaus as the genetic diversity continues to increase, meaning that after a certain genetic distance all euclidean distances in the t-SNE embedding are very similar. 
This “plateau” usually point to an embedding that gives little importance to the exact distance between clusters and instead focuses on cluster shape and spread. These statistics show that t-SNE tends to preserve global structure, shown through the high pearson coefficient, while also revealing local patterns. 
It is quite easy to distinguish viruses of different and same clades given a euclidean distance. 
In the between vs within clade violin plot for t-SNE (Figure 6C), the density of the within clade violin plot is concentrated at lower euclidean distances than the density of the between clade violin plot is, which reveals a well clustered and compartmentalized embedding of these data points. 
To quantitatively reveal this pattern, the permutation test run on this data gave a P-value of ___________________

UMAP: The embedding did translate the data differently than the rendering of the tree, where genetically similar strains were very densely packed and genetically different clades were incredibly far apart in UMAP euclidean space (Figure 4E). There were visually identifiable clusters, but the distance disparities made it hard to separate clades 9, 7, 6, and 3 from each other due to how close the clusters were.
The pattern can be seen and explained in the Hamming and genetic distance scatterplot for UMAP (Figure 5D). 
In UMAP’s scatterplot, there are two clusters of data points on the LOESS line, instead of the equal spread seen in the other 3 embeddings. 
These clusters point to a graph that places genetically similar strains incredibly close together - around 0 to 20 - but places genetically different strains almost 50 to 80 apart. 
The LOESS line reflects this increase through its exponential-like growth, and the density of these two clusters can be quantitatively shown through its Pearson Coefficient of .787. 
Euclidean distance is not a very strong measure to distinguish viruses from the same or different clades. 
In the UMAP violin plot (Figure 6D), the euclidean distance density distribution for between clusters is equal at Euclidean distances of 10 and 70. 
While the density distribution for within clusters for UMAP is concentrated almost entirely at 10, if given a euclidean distance of 10, it would be difficult to classify that strain as within the clade or between the clade it’s being compared to. 
The permutation test corroborates this, __________________


## Figure Four
:::{.out .html}
<iframe src="https://blab.github.io/cartography/FullLinkedChartClickableZika.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](FullLinkedChartClickableZika.pdf)
:::

## Figure Five

![](FullScatterplotZika.png)


## Figure Six

![](FullViolinPlotZika.png)


## SUMMARY OF RESULTS ACROSS VIRUSES

Overall, the best recapitulation of the phylogenetic clades of the four analyzed was that of MDS and UMAP, because they preserved the most local and global structure. 
However, t-SNE separated clusters much better than MDS as lots of the points were layered on top of each other in the MDS embedding (clade A1b/131K in h3n2 Influenza was impossible to see in MDS clusters 1 and 2). 

Of the 3 embeddings that reduced the Hamming distance matrix, as the focus of the algorithms shifted more towards preserving local structure over global structure, the closer the Pearson Coefficient got to 0 (MDS > UMAP > t-SNE). 
Pearson Coefficient studies the effectiveness of an embedding at preserving a relationship between genetic and euclidean distance, so for t-SNE, an algorithm that focuses primarily on exaggerating distances and clusters locally to convey patterns, the pearson coefficient is going to be closer to 0, as the data points will not adhere to a best fit line. 
For MDS, however, the embedding relies almost entirely on creating an exact 1:1 genetic:euclidean relationship, so the pearson coefficient was much higher. 
In the same vein, as the algorithms shifted towards retaining local patterns over global patterns, the disparity between the densities of the between vs within violin plots became more pronounced. 
Because the violin plots assess an embedding’s ability to distinguish between clades (how clustered the embedding is), the more exaggerated the differences between euclidean and genetic distance, the more disparate the densities are. 
t-SNE’s within clade violin plot had the most concentrated density at around 5, and its between clade violin plot had the most concentrated density at around 45. 
By comparison, the genetic distance within:between was 45:60.  

##Discussion







## Supplementary Figures and Analysis


### UMAP on leading 10 PCs:

We applied UMAP to the leading 10 PCs from the PCA plots defined below. We got qualitatively similar results to PCA, which we hypothesize to be because most of the variance in the PCA plots is explained in the first 2 PCs. Therefore, UMAP is not increasing the variance explained in the data by a significant enough amount to change the embedding qualitatively.  
##### Flu

:::{.out .html}
<iframe src="https://blab.github.io/cartography/UmapPCLinkedCladeFlu.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](UmapPCLinkedCladeFlu.pdf)
:::

##### Zika

:::{.out .html}
<iframe src="https://blab.github.io/cartography/UmapPCLinkedCladeZika.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](UmapPCLinkedCladeZika.pdf)
:::

These violin plots (the first violin plot PCA and the latter UMAP on leading 10 PCs) reveal the similar structure and density of the points found in both approaches to dimensionality reduction. While implementing UMAP on the leading 10 PCs, UMAP increases the amount of distance between clusters and concentrates datapoints within the same cluster, so the UMAP violin plot looks like better output, however, the shape of the violin plot does not significantly change, which is the qualitative measure that UMAP is not doing any clustering PCA is not already doing.  

##### Flu

![](PCAViolinPlotFlu.png)
![](UMAPPCCladeViolinPlotFlu.png)

##### Zika

![](UMAPPCCladeViolinPlotZika.png)
![](PCAViolinPlotZika.png)


### Explained Variance Plots for PCA

##### Flu 
![](ExplainedVarianceFlu.png)

##### Zika
![](ExplainedVarianceZika.png)

### PCA Full Plots

##### Flu
:::{.out .html}
<iframe src="https://blab.github.io/cartography/PCAFluBrush.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](PCAFluBrush.pdf)
:::

##### Zika
:::{.out .html}
<iframe src="https://blab.github.io/cartography/PCAZikaBrush.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](PCAZikaBrush.pdf)
:::

### TSNE Perplexity Plots:

##### Flu
:::{.out .html}
<iframe src="https://blab.github.io/cartography/TSNEPerplexityFlu.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](TSNEPerplexityFlu.pdf)
:::

##### Zika
:::{.out .html}
<iframe src="https://blab.github.io/cartography/TSNEPerplexityZika.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](TSNEPerplexityZika.pdf)
:::

### TSNE Learning Rates Plots:

##### Flu

:::{.out .html}
<iframe src="https://blab.github.io/cartography/TSNELearningRatesFlu.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](TSNELearningRatesFlu.pdf)
:::

##### Zika

:::{.out .html}
<iframe src="https://blab.github.io/cartography/TSNELearningRatesZika.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](TSNELearningRatesZika.pdf)
:::

### UMAP Nearest Neighbor and Minimum Distance Plots:

##### Flu

:::{.out .html}
<iframe src="https://blab.github.io/cartography/UMAP_neighbors_mindist_chart_Flu.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](UMAP_neighbors_mindist_chart_Flu.pdf)
:::

##### Zika

:::{.out .html}
<iframe src="https://blab.github.io/cartography/UMAP_neighbors_mindist_chart_Zika.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](UMAP_neighbors_mindist_chart_Zika.pdf)
:::

### MDS Full Plot: 

##### Flu

:::{.out .html}
<iframe src="https://blab.github.io/cartography/MDSFluBrush.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](MDSFluBrush.pdf)
:::

##### Zika

:::{.out .html}
<iframe src="https://blab.github.io/cartography/MDSZikaBrush.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](MDSZikaBrush.pdf)
:::




# Works Cited
