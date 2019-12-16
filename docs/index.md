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

#Flu 

## Description of the Problem 
<p>

Nextstrain.com can only handle diseases that have clear phylogenies, or where the virus’ ancestry is easy to track through time. However, multiple diseases, such as enterovirus, bacteria, and malaria, recombine and are therefore their phylogenies are usually split at the chromosomal level, making it hard to hypothesize the entire diseases ancestry and location throughout time. The project is to create an interactive tool that uses different embeddings of genetic sequence data via PCA, MDS, t-SNE and U-MAP to identify patterns and relationships that aren't seen in phylogenies due to this recombination issue. The overall goal would be to integrate this tool into Nextstrain in the future.

<br><br>*Recombination: Recombination is a process usually during meiosis where pieces of DNA are broken and recombined into the chromosome to produce new combinations of alleles. This process creates genetic diversity in the DNA sequences of different organisms. </p>


# METHODS:
<b>Modules necessary to run this code:</b>
<ul>
<li>Python >= 3.7</li>
<li>Docker</li>
<li>Jupyter Lab</li>
<li>BioPython</li>
<li>Pandas</li>
<li>Numpy</li>
<li>Altair</li>
<li>Seaborn</li>
<li>Scikit-Learn</li>
<li>UMAP</li>
<li>Json</li>
<li>nextstrain-augur</li>

</ul>

<br>
#Description of Methods
<p> A tree build was run using nextstrain-cli for influenza H3N2 to obtain metadata and genomic data. An aligned FASTA file was created in the results folder of a tree build from fauna <a href="https://github.com/nextstrain/fauna/tree/master/builds">here</a>. This data is available <a href="https://github.com/blab/cartography/tree/master/notebooks/Data">here</a> for free use. This data was then read and translated into 2 arrays, with the identifiers saved in strains and the sequences saved in genomes using the <a href="https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html">BioPython package SeqIO </a>.
</p>
<p>
2 different methods of transforming the data were used:
</p>
<ol>
<li> <b> Scaling and centering the data</b> </li>
<ul>
<li>Each site was treated as a different dimension, and the probability of having a certain nucleotide in that site given the frequency of that nucleotide acid at that site was reduced. This approach was modeled after the equation below where C is the dimensional matrix, M is the mean, and p is the frequency of a nucleotide at that given site </li>
<br>
<img src="C:\Users\srava\OneDrive\Pictures\Camera Roll\Equation.png" alt="Equation">
<li>Matrix M was plotted using PCA</li>
</ul>
<li> <b> Hamming distance was used to create a similarity matrix</b> </li>
</ol>
<p>
Each genome was split into separate nucleotides and compared with other nucleotides in the same site on other genomes using Hamming Distance. Only a difference between the main nucleotide pairs of A,G,C,and T was counted; Gaps (N) and other characters were not, as counting gaps would artificially cause shorter genomes to always be incredibly different compared to longer genomes (a shorter strain doesn’t necessarily mean complete genetic dissimilarity -- which is what counting gaps implies). Comparing each genome with each other genome and clustering based on their pairwise distance takes the overall structure of the multidimensional data and groups together genomes that have similar differences, which effectively clumps the data by genetic diversity. 
<br> This data was read out to a .csv file to cut processing time, and metadata was read into and merged with the genomic data (using Pandas Dataframes). The strain metadata was provided by the tree build for the influenza build (matadata_h3n2_ha), and the clade membership metadata was provided by a .json build of the influenza H3N2 tree. The similarity distance matrix was plotted using Altair (https://altair-viz.github.io/) and colored by clade assignment. 
<br>
The similarity matrix was transformed for plotting by using 4 different dimensionality reduction techniques ordered by algorithmic complexity: 
<ul> 
<li> PCA (scikit-Learn)
<li> MDS (scikit-Learn)
<li> t-SNE (scikit-Learn)
<li> UMAP (umap-Learn)
</ul>
This data was then made into a distance matrix, with each genome’s distance from another in the reduced space was plotted against the pairwise distance between the two genomes, and by using the numpy method triu_indices, which only stores values of the matrix’s upper triangle (as distance matrices are symmetrical on the diagonal), the amount of data was reduced.
<br>
Each reduction was then compared with each other in 2 different ways: 
<ol>
<li> Pairwise vs euclidean distance scatterplot with the Pearson's Coefficient as the common metric
<ul>
<li>This reveals patterns in reductions including local versus global recapitulation, effectiveness of technique at preserving accurate distances, and how clustered the data appears
</ul>
<li> Between vs Within clade distance boxplots, with the ratio of between clade distance and within clade distance used as the common metric
<ul>
<li>This reveals patterns concerning how spaced out the clusters are from each other, how well the technique clusters data accurately, and how each technique deals with outliers.
</ul>
</ol>


</p>

# RESULTS:

<p> <b> 1) Scaling and Centering the Data </b>
<br> Because PCA (Principal Component Analysis) is made to reduce multidimensional data, and not distance matrices, PCA was used to analyze the data in this format. While the ratio given by the within vs between boxplot was 47:1, which reveals a very accurately clustered set of data, the data was not truly transformed to show any new patters or information - it looked almost identical to the .json rendering of the tree. 
<iframe src="https://blab.github.io/cartography/PCAFluBrush.html" style="width: 1700px; height: 300px;" frameBorder="0"></iframe>
<iframe src="https://blab.github.io/cartography/PCABoxplot.html" style="width: 600px; height: 300px;" frameBorder="0"></iframe>
<br> <b> 2) Pairwise distance between genomes </b>
<br> The pairwise distance matrix was reduced using MDS, t-SNE, and UMAP. Comparing the pairwise vs euclidean scatterplots of these three embeddings, it can be seen that MDS's Cluster 1 and 2 Pearson Coefficient was the highest of them, sugesting the data is best correlated using MDS. The points on the scatterplot start out fairly concentrated at the origin, and spread out from there, which is an indicator that MDS is a local, rather than global, reduction technique. Locally preserving techniques tend to focus on retaining distances and structure on a smaller scale, usually within a cluster, and focuses less on retaining a "global" appearance of the data. The same pattern is observed with t-SNE to a much larger extent, as the points begin to spread out from the line of best fit much earlier on than MDS. UMAP, however, is a techinique aimed at preserving both local and global structure, and this is seen in UMAPs scatter plot. While the points do spread out as the pairwise distance gets larger, the points stay fairly concentrated in the same area. 
<iframe src="https://blab.github.io/cartography/FullScatterplot.html" style="width: 1700px; height: 300px;" frameBorder="0"></iframe>
<br> Comparing the within vs between clade boxplots for MDS, t-SNE, and UMAP, UMAP has the largest ratio of within vs between clade distance (1:7). This means that UMAP is recapitulating genetic diversity through euclidean distance with more distance than the original hamming distance matrix (which has a ratio of 1:3). MDS and t-SNE both had a ratio of 1:2. 
<iframe src="https://blab.github.io/cartography/FullBoxplot.html" style="width: 1700px; height: 300px;" frameBorder="0"></iframe>



</p>

# Tree Linked Flu Chart with PCA, MDS, TSNE, and UMAP

<iframe src="https://blab.github.io/cartography/FullLinkedChartClickable.html" style="width: 1700px; height: 850px;" frameBorder="0"></iframe>
