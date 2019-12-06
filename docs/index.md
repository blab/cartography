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

# Tree Linked Flu Chart with PCA, MDS, TSNE, and UMAP

#Flu 

## Description of the Problem 
<p>

Nextstrain.com can only handle diseases that have clear phylogenies, or where the virus’ ancestry is easy to track through time. However, multiple diseases, such as enterovirus, bacteria, and malaria, recombine and are therefore their phylogenies are usually split at the chromosomal level, making it hard to hypothesize the entire diseases ancestry and location throughout time. The project is to create an interactive tool that uses different embeddings of genetic sequence data via PCA, MDS, t-SNE and U-MAP to identify patterns and relationships that aren't seen in phylogenies due to this recombination issue. The overall goal would be to integrate this tool into Nextstrain in the future.

<br>*Recombination: Recombination is a process usually during meiosis where pieces of DNA are broken and recombined into the chromosome to produce new combinations of alleles. This process creates genetic diversity in the DNA sequences of different organisms. </p>


# RESULTS:

## Qualitative:

<p> <b> For 1) Scaling and Centering the Data </b>
<br> There didn’t seem to be any noticeable/useful patterns found in any of the embeddings
<br> <b> For 2) Pairwise distance between genomes </b>

<br> Pairwise distance matrix of genome values worked best for creating plots (show both types of plot and explain why one is better using actual numbers)
<br> By comparing each genome with each other genome and clustering based on their pairwise distance takes the overall structure of the multidimensional data and groups together genomes that have similar differences -- which means it’s clumped by either region or genetic diversity. 
<br> The less outliers there were in the dataset (pairwise distance matrix), the less accurate the recapitulation of the phylogenetic clades for the UMAP projection. 
<br> The best recapitulation of the phylogenetic clades of the four analyzed was that of MDS and t-SNE
<br> MDS had better clustering by clade and by similar regions; however, lots of the points were layered on top of each other (clade 3c3.A was impossible to see in clusters 1 and 2 and clusters 3 and 4) (as the algorithm largely makes its mapping with respect to between-object distances - local structure is preserved but not global)
t-SNE clustered genetically similar clades well, but did not place similar clusters next to each other (t-SNE has a large focus on local structure while still maintaining some global structure (presence of multiple scales of clusters))
</p>
 ## Quantitative

<p>

<br> <b> For 2) Pairwise distance between genomes </b>
	<br> The pearson coefficient for PCA for the pairwise distance vs euclidean distance for Principal Component 1 and Principal Component 2 was 0.9171.
(slope=29.26665, Y intercept=-209.08323, rvalue=0.91707, pvalue=0.0)

	<br> The pearson coefficient for MDS for the pairwise distance vs euclidean distance for clusters 1 and 2 was .8768.
(slope=0.52104, Y intercept=-2.00793, rvalue=0.87680, pvalue=0.0)

<br> The pearson coefficient for t-SNE for the pairwise distance vs euclidean distance for clusters 1 and 2 was .4250.
(slope=0.56694, Y intercept=37.01944, rvalue=0.42497, pvalue=0.0)

<br> The pearson coefficient for UMAP for the pairwise distance vs euclidean distance was .3247. 
slope=0.12132, Y intercept=11.38696, rvalue=0.32473, pvalue=0.0)

<br> <b> Analysis: </b>
<br> As the completixty of the algortihm got higher, the closer the Pearson Coefficient got to 0 (PCA > MDS > t-SNE > UMAP)
<br> As the complexity of the algorithm increased, the slove of the linear regression line (where pairwise distance is on the x and Euclidean distance is on the y) also got smaller (PCA > MDS > t-SNE > UMAP)
<br> PCA and MDS’s pearson coefficients are 2X that of t-SNE and UMAP - they may recapitulate the overall structure of the data more accurately than t-SNE and UMAP do 
<br> Euclidean distance between samples in t-SNE embeddings correlate with genetic distance with a Pearson’s R of 0.4 (p < 0.001) while UMAP embeddings only correlate with a Pearson’s R of 0.3
</p>

<iframe src="https://blab.github.io/cartography/FullLinkedChartClickable.html" style="width: 1700px; height: 850px;" frameBorder="0"></iframe>
