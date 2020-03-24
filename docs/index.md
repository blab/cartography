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

Nextstrain.com can only handle diseases that have clear phylogenies, or where the virus’ ancestry is easy to track through time. However, multiple diseases, such as enterovirus, bacteria, and malaria, recombine and are therefore their phylogenies are usually split at the chromosomal level, making it hard to hypothesize the entire diseases ancestry and location throughout time. The project is to create an interactive tool that uses different embeddings of genetic sequence data via PCA, MDS, t-SNE and U-MAP to identify patterns and relationships that aren't seen in phylogenies due to this recombination issue. The overall goal would be to integrate this tool into Nextstrain in the future.

Recombination: Recombination is a process usually during meiosis where pieces of DNA are broken and recombined into the chromosome to produce new combinations of alleles. This process creates genetic diversity in the DNA sequences of different organisms.

** PCA, MDS, t-SNE, U-MAP (Should I define each other these/ link their documentation?)

# METHODS:

### Materials:
The Modules and dependencies necessary to run this code well as instructions are on Cartography’s github (https://github.com/blab/cartography), where cartography.yml already has the modules below to install.
- pip
- seaborn
- altair
- pandas
- numpy
- biopython
- python>=3.6
- scikit-learn
- umap-learn
- jsonschema
- jupyterlab
- nextstrain-augur

[Docker](https://www.docker.com/) is also needed to run the tree builds, but to install this a computer must have certain software. More detailed instructions can be found on the website linked above. To run the code on Cartography, however, Docker is not necessary. 

### Methods:
A tree build was run using nextstrain-cli for influenza H3N2 to obtain genomic data. An aligned FASTA file was created in the results folder of the flu tree build. The build can be found at [Fauna](https://github.com/nextstrain/fauna/tree/master/builds) under "FLU" (instructions on running the build are also housed there). The created FASTA file is available [here](https://github.com/blab/cartography/tree/master/notebooks/Data) for free use. This file was then read and translated into 2 arrays using the [BioPython package SeqIO](https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html), where one array held the strain identifiers, and the other the full genome corresponding to the strain identifier.

Two different methods of transforming the data were used; Scaling and centering the data, and a Hamming distance similarity matrix. For Scaling and Centering the data, each site on the genome was treated as a different dimension, and the probability (called p) of having a certain nucleotide in that site given the frequency of that nucleotide acid at that site among all genomes was plotted. This approach was modeled after the equation below where C is the dimensional matrix, M is the mean, and p is the frequency of a nucleotide at that given site:

![](PaperImages/ScalingAndCenteringEquation.png)

Matrix M was plotted using [scikit-learn’s PCA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html) .

The second approach used Hamming distance to create a similarity matrix. Each genome was split into separate nucleotides and compared with other nucleotides in the same site on other genomes. Only a difference between the main nucleotide pairs (AGCT) was counted -- gaps (N) were not. This is because some sequences were significantly shorter than others, and a shorter strain does not necessarily mean complete genetic dissimilarity, which is what counting gaps implied. An illustration follows below to demonstrate this process. 

![](PaperImages/SimilarityMatrixExplanation.png)

The similarity matrix was read out to a .csv file to cut processing time. The similarity distance matrix was reduced through MDS, t-SNE, and UMAP, plotted using [Altair](https://altair-viz.github.io/) ,and colored by clade assignment. Clade membership metadata was provided by a .json build of the influenza H3N2 tree (the build can be found at https://github.com/blab/cartography/tree/master/notebooks/Data) . The 3 different dimensionality reduction techniques are ordered below by algorithmic complexity: 
- MDS (scikit-Learn)
- t-SNE (scikit-Learn)
- UMAP (umap-Learn)

To further analyze the embeddings’ ability to accurately capture the multidimensional data, two separate plots were made: pairwise vs euclidean distance scatterplots with a LOESS best fit line, and within vs between clade violin plots per embedding.

Pairwise vs euclidean distance scatterplots:
 
The similarity matrix’s upper triangle values are made into a flattened numpy array using [numpy.triu_indices()](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.triu_indices.html) and joined with a flattened numpy array of euclidean distances between each embedding’s coordinate points in a Pandas Dataframe. This dataframe is plotted using seaborn and [statsplots’ LOESS function](https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html) . This data was made into a distance matrix as well, where each genome’s distance from another in the reduced space was plotted against the pairwise distance between the two genomes. The pairwise distance was on the x axis, and the euclidean distance was on the y axis. Linear regression data was calculated from the Pandas Dataframe using [linregress](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html).

Between vs Within clade Violin plots:

The matrix of euclidean distances for each embedding was flattened, and each comparison was labeled as a “within clade” or “between clade” comparison using the clade assignments from the .json build of the tree. Violin plots were made using [seaborn](https://seaborn.pydata.org/) , separated by clade status and euclidean distance on the y axis.  

# RESULTS:

### Qualitative:

Scaling and Centering the Data

Influenza:
- Because PCA (Principal Component Analysis) reduces multidimensional data and not distance matrices, PCA was used to analyze the data in this format. While the ratio given by the within vs between violin plot was 24:1 and a positive R^2 correlation of .57, revealing a tightly clustered set of data, the data was not transformed to show any new pattern or information, and clustered the data almost identically to the .json rendering of the tree.
<iframe src="C:/Users/srava/BedfordProjects/cartography/docs/PCAFluBrush.html" style="width: 800px; height: 400px;" frameBorder="0"></iframe>
![](PCAViolinPlotFlu.png)
![](PCAScatterplotLOESSFlu.png)
Zika:
- The Ratio given by the within vs between violin plot was 1:1, revealing data not clustered at all. The chart reflects this, having no pattern whatsoever. 
Pairwise distance between genomes

Pairwise distance matrix of genome values worked best for creating plots (show both types of plot and explain why one is better using actual numbers)
By comparing each genome with each other genome and clustering based on their pairwise distance takes the overall structure of the multidimensional data and groups together genomes that have similar differences -- which means it’s clumped by either region or genetic diversity. 
The less outliers there were in the dataset (pairwise distance matrix), the less accurate the recapitulation of the phylogenetic clades for the UMAP projection. 
The best recapitulation of the phylogenetic clades of the four analyzed was that of MDS and t-SNE
MDS had better clustering by clade and by similar regions; however, lots of the points were layered on top of each other (clade 3c3.A was impossible to see in clusters 1 and 2 and clusters 3 and 4) (as the algorithm largely makes its mapping with respect to between-object distances - local structure is preserved but not global)
t-SNE clustered genetically similar clades well, but did not place similar clusters next to each other (t-SNE has a large focus on local structure while still maintaining some global structure (presence of multiple scales of clusters))
Quantitative

For 2) Pairwise distance between genomes
	The pearson coefficient for PCA for the pairwise distance vs euclidean distance for Principal Component 1 and Principal Component 2 was 0.9171.
(slope=29.26665, Y intercept=-209.08323, rvalue=0.91707, pvalue=0.0)

	The pearson coefficient for MDS for the pairwise distance vs euclidean distance for clusters 1 and 2 was .8768.
(slope=0.52104, Y intercept=-2.00793, rvalue=0.87680, pvalue=0.0)

The pearson coefficient for t-SNE for the pairwise distance vs euclidean distance for clusters 1 and 2 was .4250.
(slope=0.56694, Y intercept=37.01944, rvalue=0.42497, pvalue=0.0)

The pearson coefficient for UMAP for the pairwise distance vs euclidean distance was .3247. 
slope=0.12132, Y intercept=11.38696, rvalue=0.32473, pvalue=0.0)

Analysis:
As the completixty of the algortihm got higher, the closer the Pearson Coefficient got to 0 (PCA > MDS > t-SNE > UMAP)
As the complexity of the algorithm increased, the slope of the linear regression line (where pairwise distance is on the x and Euclidean distance is on the y) also got smaller (PCA > MDS > t-SNE > UMAP)
PCA and MDS’s pearson coefficients are 2X that of t-SNE and UMAP - they may recapitulate the overall structure of the data more accurately than t-SNE and UMAP do 
 Euclidean distance between samples in t-SNE embeddings correlate with genetic distance with a Pearson’s R of 0.4 (p < 0.001) while UMAP embeddings only correlate with a Pearson’s R of 0.3
The 









#Flu 

## Description of the Problem 
<p>

Most diseases can be modeled in a transmission tree, a modeling approach that tracks mutations in disease samples back to a single common ancestor. The ground assumption for a transmission tree is that the mutations present in the samples can be aligned and standardized - this creates a clear phylogeny. Nextstrain.com can only handle diseases that have clear phylogenies, or where the virus’ ancestry is easy to track through time. However, recombinant diseases such as enterovirus, bacteria, and malaria switch chromosome pairs when they combine, making it impossible to quantify a genetic sequence’s mutations in reference to another. To combat this problem, most epidemiologist split their phylogenies at the chromosomal level, but this makes it incredibly difficult to hypothesize the entire diseases ancestry and location throughout time. The project is to create an interactive tool that uses different embeddings of genetic sequence data via PCA, MDS, t-SNE and U-MAP to identify patterns and relationships that aren't seen in phylogenies due to this recombination issue. The overall goal would be to integrate this tool into Nextstrain in the future.

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
<li> PCA (scikit-Learn)</li>
<li> MDS (scikit-Learn)</li>
<li> t-SNE (scikit-Learn)</li>
<li> UMAP (umap-Learn)</li>
</ul>
This data was then made into a distance matrix, with each genome’s distance from another in the reduced space was plotted against the pairwise distance between the two genomes, and by using the numpy method triu_indices, which only stores values of the matrix’s upper triangle (as distance matrices are symmetrical on the diagonal), the amount of data was reduced.
<br>
Each reduction was then compared with each other in 2 different ways: 
<ol>
<li> Pairwise vs euclidean distance scatterplot with the Pearson's Coefficient as the common metric</li>
<ul>
<li>This reveals patterns in reductions including local versus global recapitulation, effectiveness of technique at preserving accurate distances, and how clustered the data appears.</li>
</ul>
<li> Between vs Within clade distance boxplots, with the ratio of between clade distance and within clade distance used as the common metric</li>
<ul>
<li>This reveals patterns concerning how spaced out the clusters are from each other, how well the technique clusters data accurately, and how each technique deals with outliers.</li>
</ul>
</ol>


</p>

# RESULTS:

<p> <b> 1) Scaling and Centering the Data </b>
<br> Because PCA (Principal Component Analysis) is made to reduce multidimensional data, and not distance matrices, PCA was used to analyze the data in this format. While the ratio given by the within vs between boxplot was 47:1, which reveals a very accurately clustered set of data, the data was not truly transformed to show any new patters or information - it looked almost identical to the .json rendering of the tree. 
<iframe src="https://blab.github.io/cartography/PCAFluBrush.html" style="width: 1700px; height: 300px;" frameBorder="0"></iframe>
<iframe src="https://blab.github.io/cartography/PCABoxplot.html" style="width: 600px; height: 300px;" frameBorder="0"></iframe>
<b> 2) Pairwise distance between genomes </b>
<br> The pairwise distance matrix was reduced using MDS, t-SNE, and UMAP. Comparing the pairwise vs euclidean scatterplots of these three embeddings (excluding PCA, as it doesn't reduce similarity matrix data), it can be seen that MDS's Cluster 1 and 2 Pearson Coefficient was the highest of them (r = .595), sugesting the data is best correlated using MDS. The points on the scatterplot start out fairly concentrated at the origin, and spread out from there, which is an indicator that MDS is a local, rather than global, reduction technique. Locally preserving techniques tend to focus on retaining distances and structure on a smaller scale, usually within a cluster, and focuses less on retaining a "global" appearance of the data. The same pattern is observed with t-SNE to a much larger extent, as the points begin to spread out from the line of best fit much earlier on than MDS. UMAP, however, is a techinique aimed at preserving both local and global structure, and this is seen in UMAPs scatter plot. While the points do spread out as the pairwise distance gets larger, the points stay fairly concentrated in the same area. 
<iframe src="https://blab.github.io/cartography/FullScatterplot.html" style="width: 1700px; height: 300px;" frameBorder="0"></iframe>
<br> Comparing the within vs between clade boxplots for MDS, t-SNE, and UMAP, UMAP has the largest ratio of within vs between clade distance (1:7). This means that UMAP is recapitulating genetic diversity through euclidean distance with more distance than the original hamming distance matrix (which has a ratio of 1:3). MDS and t-SNE both had a ratio of 1:2. 
<iframe src="https://blab.github.io/cartography/FullBoxplot.html" style="width: 1700px; height: 400px;" frameBorder="0"></iframe>



</p>

# Tree Linked Flu Chart with PCA, MDS, TSNE, and UMAP

<iframe src="https://blab.github.io/cartography/FullLinkedChartClickable.html" style="width: 1700px; height: 850px;" frameBorder="0"></iframe>
