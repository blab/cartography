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
MDS (scikit-Learn)
t-SNE (scikit-Learn)
UMAP (umap-Learn)
To further analyze the embeddings’ ability to accurately capture the multidimensional data, two separate plots were made: pairwise vs euclidean distance scatterplots with a LOESS best fit line, and within vs between clade violin plots per embedding.

Pairwise vs euclidean distance scatterplots: 
The similarity matrix’s upper triangle values are made into a flattened numpy array using [numpy.triu_indices()](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.triu_indices.html) and joined with a flattened numpy array of euclidean distances between each embedding’s coordinate points in a Pandas Dataframe. This dataframe is plotted using seaborn and [statsplots’ LOESS function](https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html) . This data was made into a distance matrix as well, where each genome’s distance from another in the reduced space was plotted against the pairwise distance between the two genomes. The pairwise distance was on the x axis, and the euclidean distance was on the y axis. Linear regression data was calculated from the Pandas Dataframe using [linregress](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html) .

Between vs Within clade Violin plots

The matrix of euclidean distances for each embedding was flattened, and each comparison was labeled as a “within clade” or “between clade” comparison using the clade assignments from the .json build of the tree. Violin plots were made using [seaborn](https://seaborn.pydata.org/) , separated by clade status and euclidean distance on the y axis.  

RESULTS:
(PCAViolinPlotFlu.png)
Qualitative:
Scaling and Centering the Data
Influenza:
- Because PCA (Principal Component Analysis) reduces multidimensional data and not distance matrices, PCA was used to analyze the data in this format. While the ratio given by the within vs between violin plot was 24:1 and a positive R^2 correlation of .57, revealing a tightly clustered set of data, the data was not transformed to show any new pattern or information, and clustered the data almost identically to the .json rendering of the tree.
<iframe src="https://blab.github.io/cartography/PCAFluBrush.html" style="width: 800px; height: 400px;" frameBorder="0"></iframe>
![]((PCAViolinPlotFlu.png)
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
