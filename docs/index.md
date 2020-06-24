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
One approach has been to split a genome into multiple phylogenies to model the evolution of the nonrecombinant fragments by using a genetic algorithm that scans strains for recombination breakpoints, quantifies and analyzes the impact of recombination at each breakpoint, and breaks the phylogeny at the most important breakpoints [@kosakovsky_pond_posada_gravenor_woelk_frost_2006].
Finding recombination breakpoints relies on detection of a recombination signal through methods such as CHIMAERA and LARD, which use split decomposition, a method which depicts parallel edges between sequences if there are conflicting phylogenetic signals in the data [@Posada13757][@10.1093/oxfordjournals.molbev.a026121][@martin_murrell_khoosal_muhire_2017].
An alternate strategy is to compare viral genomes with alternative methods that do not make the same strong assumptions as phylogenetic inference (e.g., PCA, MDS, etc.)
PCA has been used to estimate and model ancestry in unrelated individuals and plot peoples genomes to reveal patterns in national origin [@alexander_novembre_lange_2009][@GenesMirrorGeography2008], and was also used to genotype major classes of structural variants - structural differences such as deletion, duplication, and novel sequence insertion - in diverse populations to map their population stratification [@sudmant_1000_gene_paper].
PCA was also used to reveal Zika’s genetic diversity and spread in the Americas by assessing clustering of multidimensional genetic data [@H.C.109348].
MDS has been applied to H3N2 sequences to inspect relationships between all gene segments, which is closely related to the subject of this paper, for the difference that Rambaut et .al. 2008 looks at between-gene diversity rather than within-gene [@rambaut_pybus_nelson_viboud_taubenberger_holmes_2008].
PCA, t-SNE, and UMAP have all been used to describe and correlate human genomes and characteristics, which has shown that even within broader population groups, family members were projected nearer to each other in the embeddings[diaz-papkovich_anderson-trocmé_ben-eghan_gravel_2019]. 
We present a novel approach to understanding relationships among viral genomes by transforming genomic data and then using dimensionality reduction methods such as PCA, MDS, t-SNE, and UMAP. 
We investigate the degree to which this method can recapitulate known phylogenetic relationships for viruses whose genomes are phylogenetically tractable (we used influenza H3N2 and Zika).
We apply this method to viruses whose genomes are known to undergo substantial recombination, such as MERS and SARS-CoV-2 to assess how well each method is able to reconstruct previously identified biologically-meaningful clusters.

Recombination: occurs when at least two viral genomes co-infect the same host cell and exchange genetic segments. 
Shuffling/reassortment, a particular type of recombination, occurs in viruses with segmented genomes, which by interchanging complete genome segments, gives rise to new segment combinations [@pérez-losada_arenas_galán_palero_gonzález-candelas_2014].


# METHODS:

### Materials:

The analysis environment can be recreated using conda and all installation instructions are available on [Cartography’s github](https://github.com/blab/cartography), where cartography.yml already has the modules below to install.

### Methods:
We analyzed Influenza A/H3N2 and Zika genomes and created a FASTA file of multiple sequence alignments with MAFFT v7.407 [@Katoh2002] via augur align [@Hadfield2018] and phylogenies with IQ-TREE v1.6.10 [@Nguyen2014] via augur tree version 5.4.1.

We used two different methods of transforming the data; Scaling and centering the data, and a Hamming distance similarity matrix.
For Scaling and Centering the data, we performed PCA on the matrix of nucleotides from the multiple sequence alignment using scikit-learn [@jolliffe_cadima_2016].

The second approach we used Hamming distance was to create a similarity matrix. 
By comparing every genome with every other genome and clustering based on their pairwise distance, distance-based methods take the overall structure of the multidimensional data and groups together genomes that have similar differences.
This means the data is clustered by genetic diversity (in a phylogenetic tree genetic diversity is categorized using clades).
Each genome was split into separate nucleotides and compared with other nucleotides in the same site on other genomes.
We only counted a difference between the main nucleotide pairs (AGCT) -- gaps (N) were not.
This is because some sequences were significantly shorter than others, and a shorter strain does not necessarily mean complete genetic dissimilarity, which is what counting gaps implied. 

We reduced the similarity distance matrix through MDS, t-SNE, and UMAP, plotted using [Altair](https://altair-viz.github.io/) ,and colored by clade assignment.
Clade membership metadata was provided by a .json build of the influenza H3N2 tree (the build can be found at https://github.com/blab/cartography/tree/master/notebooks/Data).
The 3 different dimensionality reduction techniques are ordered below by publication date: 
- [MDS](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.MDS.html)
- [t-SNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)
- [UMAP](https://umap-learn.readthedocs.io/en/latest/)

To further analyze the embeddings’ ability to accurately capture the multidimensional data, we made two separate plots: hamming vs euclidean distance scatterplots with a LOESS best fit line, and within vs between clade violin plots per embedding.

Hamming distance vs euclidean distance scatterplots:
 
Hamming distance vs Euclidean distance plots assess the local and global structure of the embedding as well as assess the overall strength of the embedding's recapitulation.
The Hamming distance of nucleotide sequences is plotted on the x axis, and the euclidean distance between the points in the embedding are plotted on the y axis.
By plotting these distance measurements, we can observe how correlated the dataset is.
The higher the correlation, the better a function can describe the relationship between the Hamming distance value and the euclidean distance value. 
In this way, constant correlation in a plot reveals that the embedding tends to capture and retain local patterns rather than global, and splayed points point to a global structure preservation over local. 
Therefore, the closer the Pearson Coefficient is to 1, the better the embedding is at preserving genetic dissimilarity in euclidean space.
The LOESS line drawn through the plot assesses the best fit function for the embedding.



Between vs Within clade Violin plots:

The Between vs Within clade Violin plots visually represent how well Euclidean distances can distinguish virus genomes from different clades, or shows the probability of a certain Euclidean distance being used to classify the within vs between status of a given pair of genomes. 
The density of the violin plot at a specific distance gives the relative probability of two strains with that distance being in same or different clades.
The median to median ratio is an indicator of how well the embedding clusters the data.
The larger the ratio between the medians of between vs within clade violin plots, the better the embedding is at clustering and compartmentalizing data into their clades.
To create this plot, the matrix of euclidean distances for each embedding was flattened, and each comparison was labeled as a “within clade” or “between clade” comparison using the clade assignments from the .json build of the tree.
Violin plots were made using [seaborn](https://seaborn.pydata.org/) , separated by clade status and euclidean distance on the y axis.  

# RESULTS:

### EXPECTATIONS FOR PCA, MDS, t-SNE, and UMAP

Principal Component Analysis (PCA) reduces multidimensional data, increasing interpretability while minimizing information loss.
PCA relies on linear assumptions, does not affect the scale of the data, and does not normalize the data as part of the algorithm. 
PCA is almost entirely focused on retaining the global structure and variance of the data, and therefore one of its limitations is revealing patterns locally. 
PCA is not an algorithm to be used on a similarity matrix, and is instead intended for transformed and normalized multidimensional data.
In the context of this paper, PCA will be used on transformed and normalized genetic data and not on the similarity matrix described above.

Multidimensional Scaling (MDS) refers to statistical techniques that reduce the complexity of a data set by quantifying similarity judgments, which increases the interpretability of local relational structures mired in the dataset.
A limitation to MDS is that only one symmetric matrix is allowed as input, and the scale of measurement is non-numerical.
MDS preserves global patterns over local, but the algorithm’s importance on translating dissimilarity to distance does preserve some larger patterns in local structure as well.
In the context of this paper, MDS will cluster data into sparse “sections” of the map while not creating actual clusters. 

t-distributed Stochastic Neighbor Embedding (t-SNE) visualizes high-dimensional data by giving each datapoint a location in a 2 to 3 dimensional map.
t-SNE is focused largely on local structure over global structure, and t-SNE’s projection of clusters and distances between clusters are not analogous to dissimilarity - in other words, t-SNE focuses heavily on projecting similarity rather than dissimilarity [@maaten2008visualizing].
Because t-SNE reduces data’s dimensionality based on local properties of data, data with intrinsically high dimensional structure will not be projected accurately.
In the context of this paper, t-SNE will create tight clusters that clearly indicate genetic similarity, but will not create an accurate global picture of the data.

UMAP (Uniform Manifold Approximation and Projection) is a manifold learning technique for dimension reduction based in Riemannian geometry and algebraic topology [@lel2018umap].
Some limitations include its lack of maturity - this novel technique does not have firmly established or robust practices and libraries to use UMAP best.
In the context of this paper, UMAP will reveal a tightly clustered set of data that retains both the global structure of the data and the clusters and similarities present at the local level. 

# EXPECTATIONS FOR VIRUSES:

## Influenza:

H3N2 Influenza in this project is used as a proof of concept as influenza only reassorts and does not recombine. 
We use H3N2's HA sequences as they have a relatively high mutation rate compared to the other gene segments, it encodes a protein that is a target of human immunity, and has traditionally been used for analysis of influenza evolution. 
As these sequences are biologically relevant, short, and do not recombine, the genomes can be reasonably assigned to phylogenetic clades. 
Therefore, it can be assumed that H3N2 HA influenza is a good test case for Cartography.

## SUMMARY OF RESULTS FOR INFLUENZA 

PCA: The data was normalized and scaled, and PCA was run on it (Figure 1).
The PCA pairwise vs euclidean distance plot gave a pearson coefficient of .739, which was the closest coefficient to 1 of the 4 embeddings analyzed.
The data points followed the LOESS line in a fairly linear fashion, which, given that PCA retains global patterns and therefore genetic distances should not differ widely from euclidean distance in the plot, upholds preexisting beliefs about the algorithm (Figure 2).
The ratio between the medians given by the within vs between violin plot for PCA was 24:1, revealing tightly clustered and compartmentalized data (Figure 3).
However, the data did not seem to show any new patterns of information and clustered the data points almost identically to the .json rendering of the tree; similar distance and placement of clade clusters corroborates this hypothesis.


MDS: MDS reduced the pairwise distance matrix (Figure 1).
The MDS pairwise vs euclidean distance plot gave a pearson coefficient of .615, which was the highest of the 3 embeddings that reduced pairwise distance matrix data.
The data points followed the LOESS line in a fairly linear fashion, which, given that MDS favors retaining global patterns over local patterns and therefore genetic distances should not differ widely from euclidean distance in the plot, upholds preexisting beliefs about the algorithm (Figure 2).
The ratio between the medians given by the within vs between violin plot for MDS was 4:1, revealing fairly tightly clustered data.
The data points from different clades did overlap each other, but given that MDS retains global structure, this result is expected.
The higher density of the MDS  between vs within violin plots at lower euclidean distances corroborates this pattern, as it reveals that it is harder to distinguish the between vs within clade measure at lower MDS euclidean distances (Figure 3).
The embedding did translate the data differently than the rendering of the tree, giving more euclidean space between points that were more genetically diverse.
This pattern is clearly shown in the pairwise vs euclidean graph for MDS, as the plateau in the lower genetic distances quickly gives way to a steep incline as euclidean and genetic distance begin to diverge increasingly.


t-SNE: t-SNE reduced the pairwise distance matrix (Figure 1).
The t-SNE pairwise vs euclidean distance plot gave a pearson coefficient of .307, which was the lowest of the 3 embeddings that reduced pairwise distance matrix data.
The data points splayed out from the LOESS line as genetic distance increased along the x axis, which, given that t-SNE favors retaining local patterns over global patterns and therefore euclidean and genetic distance begin to diverge increasingly, upholds preexisting beliefs about the algorithm (Figure 2).
The ratio between the medians given by the within vs between violin plot for t-SNE was 3:1, revealing fairly clustered data.
The data points from different clades did not overlap each other as t-SNE is a locally focused embedding.
The density of the within clade violin plot is at lower euclidean distances than the density of the between clade violin plot, revealing that it is quite easy to distinguish the between vs within clade measure given a euclidean distance (this overall reveals a well clustered and compartmentalized embedding of these data points) (Figure 3).
The embedding did translate the data differently than the rendering of the tree, giving more euclidean space between points that were more genetically diverse.
The t-SNE embedding translated the data and found local clusters previously unseen; one of the clades (A1b/131K) was split into three clusters that, when looking at NextStrain’s updated Ha phylogeny, slices the clade in three clean pieces. Looking at the nextstrain clade, this “slicing” in the dimensionality reduction revealed patterns in the tree that were not seen, and this foreshadows Cartography’s possible competency with pinpointing future clade break points.


UMAP: UMAP reduced the pairwise distance matrix (Figure 1).
The UMAP pairwise vs euclidean distance plot gave a pearson coefficient of .471.
The data points did splay out from the LOESS line, but far less so compared to t-SNE - the more distinguishing characteristics were the clusters shown in the scatter plot itself, revealing UMAP’s approach to visualizing local and global patterns equally (Figure 2).
The ratio between the medians given by the within vs between violin plot for UMAP was 9:1, revealing tightly clustered data.
The data points from different clades were far apart and incredibly densely packed.
The density of the within clade violin plot is at much lower euclidean distances than the density of the between clade violin plot, revealing that it is extremely easy to discern the between vs within clade measure given a euclidean distance (corroborating the observation of UMAP’s densely packed and compartmentalized embedding of these data points) (Figure 3).
The embedding did translate the data differently than the rendering of the tree, giving much more euclidean space between points that were more genetically diverse and packing genetically similar strains very close together, if not stacked.
This pattern is clearly shown in the pairwise vs euclidean graph for UMAP, as the LOESS line begins to plateau as the genetic distance gets larger and larger, revealing that UMAP sets a threshold for genetic diversity and keeps points above that threshold very far apart in euclidean space.
This is different from other embeddings, because PCA, MDS, and t-SNE do not set a threshold for genetic diversity, and instead scale the data based on the largest genetic distance disparity (this makes UMAP’s embedding so densely clustered).


## Figure One
:::{.out .html}
<iframe src="https://blab.github.io/cartography/FullLinkedChartClickableFlu.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](docs/FullLinkedChartClickableFlu.pdf)
:::


## Figure Two
:::{.out .html}
![](FullScatterplotFlu.png)
:::
:::{.out .pdf}
![](docs/FullScatterplotFlu.png)
:::
## Figure Three
:::{.out .html}
![](FullViolinPlotFlu.png)
:::
:::{.out .pdf}
![](docs/FullViolinPlotFlu.png)
:::


## Zika:

Zika: Zika in this project is used as a test case. 
While H3N2 Influenza is a globally distributed virus that has caused infections seasonally for decades, Zika is a fairly new human pathogenic virus that has a restricted geographic distribution that recapitulates the patterns of viral transmission. 
Therefore, Zika is better compartamentalized by region than by clade. 
It can be reasonably assumed that zika is a good test case for Cartography.

## SUMMARY OF RESULTS FOR ZIKA: 

PCA: The data was normalized and scaled, and PCA was run on it (Figure 4). 
The PCA pairwise vs euclidean distance plot gave a pearson coefficient of .005, which was the lowest coefficient of the 4 embeddings analyzed. 
The data points did not have any correlation to each other, meaning the data showed no signs of clustering at all (Figure 5). 
The ratio between the medians given by the within vs between violin plot for PCA was 1:1, meaning there is no distinction between between clade and within clade distance, corroborating previous conclusions (Figure 6).

MDS: MDS reduced the pairwise distance matrix (Figure 4). 
The MDS pairwise vs euclidean distance plot gave a pearson coefficient of .878, which was the highest of the 4 embeddings. 
The data points followed the LOESS line in an extremely linear fashion, which, given that MDS favors retaining global patterns over local patterns and therefore genetic distances should not differ widely from euclidean distance in the plot, upholds preexisting beliefs about the algorithm (Figure 5). 
The ratio between the medians given by the within vs between violin plot for MDS was 7:1, revealing very tightly clustered data. 
Some data points from different clades did overlap each other, but given that MDS retains global structure, this result is expected.
The density of the within clade violin plot is at lower euclidean distances than the density of the between clade violin plot, revealing that it is quite easy to distinguish the between vs within clade measure given a euclidean distance (Figure 6). 
There strong correlation between in the axes in the pairwise vs euclidean graph for MDS reveal an embedding that translates hamming distance to euclidean distance without much variability. 

t-SNE: t-SNE reduced the pairwise distance matrix (Figure 4). 
The t-SNE pairwise vs euclidean distance plot gave a pearson coefficient of .530, which was the lowest of the 3 embeddings that reduced pairwise distance matrix data. 
The data points splayed out from the LOESS line as genetic distance increased along the x axis, which, given that t-SNE favors retaining local patterns over global patterns and therefore euclidean and genetic distance begin to diverge increasingly, upholds preexisting beliefs about the algorithm (Figure 5). 
The Euclidean vs Genetic distance scatter plot grows in a linear fashion until it plateaus, meaning that after a certain genetic distance all euclidean distances are similar. 
The ratio between the medians given by the within vs between violin plot for t-SNE was 7:1, revealing very clustered data. 
The data points from different clades did not overlap each other very much as t-SNE is a locally focused embedding. 
The density of the within clade violin plot is at lower euclidean distances than the density of the between clade violin plot, revealing that it is quite easy to distinguish the between vs within clade measure given a euclidean distance (this overall reveals a well clustered and compartmentalized embedding of these data points) (Figure 6). 
The embedding did translate the data differently than the rendering of the tree, giving less euclidean space between points that were more genetically diverse, and creating a “threshold” for how much genetic diversity would impact euclidean distance. 

UMAP: UMAP reduced the pairwise distance matrix (Figure 4). 
The UMAP pairwise vs euclidean distance plot gave a pearson coefficient of .753. 
The data points did splay out from the LOESS line, but far less so compared to t-SNE - the more distinguishing characteristics were the clusters shown in the scatter plot itself that highlights UMAP’s stepwise approach to clustering - UMAP places much more distance between genetically diverse strains than the global embedding as UMAP visualizes local and global patterns equally (Figure 5). 
The ratio between the medians given by the within vs between violin plot for UMAP was 7:1, revealing tightly clustered data. 
The data points from different clades were far apart and incredibly densely packed when similar. The density of the within clade violin plot is at much lower euclidean distances than the density of the between clade violin plot, revealing that it is quite easy to discern the between vs within clade measure given a euclidean distance (Figure 6).
The embedding did translate the data differently than the rendering of the tree, giving much more euclidean space between points that were more genetically diverse and packing genetically similar strains very close together, if not stacked. 
This is different from other embeddings, because PCA, MDS, and t-SNE scale the data based on the largest genetic distance disparity instead of UMAP’s ability to scale more euclidean distance between points very genetically dissimilar (scaling euclidean distance exponentially compared to genetic distance).

## Figure Four
:::{.out .html}
<iframe src="https://blab.github.io/cartography/FullLinkedChartClickableZika.html" style="width: 1200px; height: 400px;" frameBorder="0"></iframe>
:::
:::{.out .pdf}
![](docs/FullLinkedChartClickableZika.pdf)
:::

## Figure Five
:::{.out .html}
![](FullScatterplotZika.png)
:::
:::{.out .pdf}
![](docs/FullScatterplotZika.png)
:::
## Figure Six
:::{.out .html}
![](FullViolinPlotZika.png)
:::
:::{.out .pdf}
![](docs/FullViolinPlotZika.png)
:::

## SUMMARY OF RESULTS ACROSS VIRUSES

Overall, the best recapitulation of the phylogenetic clades of the four analyzed was that of MDS and UMAP, because they preserved the most local and global structure. 
However, t-SNE separated clusters much better than MDS as lots of the points were layered on top of each other in the MDS embedding (clade A1b/131K in H3N2 Influenza was impossible to see in MDS clusters 1 and 2). 

Of the 3 embeddings that reduced the pairwise distance matrix, as the focus of the algorithms shifted more towards preserving local structure over global structure, the closer the Pearson Coefficient got to 0 (MDS > UMAP > t-SNE). 
Pearson Coefficient studies the effectiveness of an embedding at preserving a relationship between genetic and euclidean distance, so for t-SNE, an algorithm that focuses primarily on exaggerating distances and clusters locally to convey patterns, the pearson coefficient is going to be closer to 0, as the data points will not adhere to a best fit line. 
For MDS, however, the embedding relies almost entirely on creating an exact 1:1 genetic:euclidean relationship, so the pearson coefficient was much higher. 
In the same vein, as the algorithms shifted towards retaining local patterns over global patterns, the disparity between the densities of the between vs within violin plots became more pronounced. 
Because the violin plots assess an embedding’s ability to distinguish between clades (how clustered the embedding is), the more exaggerated the differences between euclidean and genetic distance, the more disparate the densities are. 
t-SNE’s within clade violin plot had the most concentrated density at around 5, and its between clade violin plot had the most concentrated density at around 45. 
By comparison, the genetic distance within:between was 45:60.  


# Works Cited
