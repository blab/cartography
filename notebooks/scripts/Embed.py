"""
takes the similarity matrix and the strain and genome CSV files and outputs the dataframe to be visualized using altair/matplotlib
"""
import pandas as pd
import altair as alt
import numpy as np
from scipy.spatial.distance import squareform, pdist
from Bio import SeqIO
import seaborn as sns
import re
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from augur.utils import json_to_tree
import json
from sklearn.manifold import MDS
from sklearn.manifold import TSNE
import umap
from scipy.stats import linregress
from pathlib import Path
import statsmodels
import statistics
import matplotlib.pyplot as plt
from Helpers import get_euclidean_data_frame, get_hamming_distances, linking_tree_with_plots_brush
from Helpers import linking_tree_with_plots_clickable
from Helpers import scatterplot_xyvalues, scatterplot_tooltips, scatterplot_with_tooltip_interactive
from selenium.webdriver import Chrome 
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from Helpers import get_y_positions
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--method", required=True, help="FASTA file of sequences")
parser.add_argument("--parameters", required=False, help="parameter options: perplexity, learning rate, nearest neighbors, min_dist")
parser.add_argument("--disease_name", required=True, help="name of disease")
parser.add_argument("--output-node-data", required=True, help="output formats: json or node js for auspice")
	
args = parser.parse_args()


similarity_matrix = pd.read_csv('../Dataframes/distance_matrix_' + args.disease_name + '.csv', index_col=0)

strains_df = pd.read_csv('../Dataframes/strains_dataframe' + args.disease_name + '.csv') 
 
genomes_df = pd.read_csv('../Dataframes/genomes_dataframe' + args.disease_name + '.csv')

node_df = pd.read_csv('../Dataframes/node_dataframe' + args.disease_name + '.csv')

genomes = genomes_df["strain"].values.tolist()

strains = strains_df["strain"].values.tolist()



def PCA():
	numbers = genomes[:]
	for i in range(0,len(genomes)):
		numbers[i] = re.sub(r'[^AGCT]', '5', numbers[i])
		numbers[i] = list(numbers[i].replace('A','1').replace('G','2').replace('C', '3').replace('T','4'))
		numbers[i] = [int(j) for j in numbers[i]]
	genomes_df = pd.DataFrame(numbers)
	genomes_df.columns = ["Site " + str(k) for k in range(0,len(numbers[i]))]

	#performing PCA on my pandas dataframe 
	pca = PCA(n_components=10,svd_solver='full') #can specify n, since with no prior knowledge, I use None
	principalComponents = pca.fit_transform(genomes_df)

	# Create a data frame from the PCA embedding.
	principalDf = pd.DataFrame(data = principalComponents, columns = ["PCA" + str(i) for i in range(1,11)])

	# Annotate rows by their original strain names. PCA rows are in the same order as
	# the `genomes` rows which are in the same order as the `strains` rows.
	principalDf["strain"] = strains

	#dataframe for explained variance plot
	df = pd.concat([pd.DataFrame(np.arange(1,11)), pd.DataFrame([round(pca.explained_variance_ratio_[i],4) for i in range(0,len(pca.explained_variance_ratio_))])], axis = 1)
	df.columns = ['principal components','explained variance']

	#merged pca dataframe
	merged_pca_df = principalDf.merge(node_df, on="strain")
	
	merged_pca_df.to_csv('../Dataframes/PCA_dataframe' + args.disease_name + '.csv', index = False)

	violin_plot_and_scatterplot_code(merged_pca_df, "PCA1", "PCA2", "PCA")
def MDS():
	embedding = MDS(n_components=10,metric=True,dissimilarity='precomputed')
	X_transformed_mds = embedding.fit_transform(similarity_matrix)
	raw_stress = embedding.stress_
	normalized_stress = np.sqrt(raw_stress /((similarity_matrix.values.ravel() ** 2).sum() / 2)
	MDS_df = pd.DataFrame(X_transformed_mds,columns=['MDS' + str(i) for i in range(1,11)])
	# Annotate rows by their original strain names. 
	MDS_df["strain"] = similarity_matrix.index
	merged_mds_df = MDS_df.merge(node_df, on="strain")
	
	merged_mds_df.to_csv('../Dataframes/MDS_dataframe' + args.disease_name + '.csv', index = False)
	
	violin_plot_and_scatterplot_code(merged_mds_df, "MDS1", "MDS2", "MDS")

def TSNE():
	embedding = TSNE(n_components=2,metric='precomputed',perplexity = 25.95)
	X_transformed_tsne = embedding.fit_transform(similarity_matrix)
	TSNE_df = pd.DataFrame(X_transformed_tsne,columns=['TSNE' + str(i) for i in range(1,3)])
	TSNE_df["strain"] = similarity_matrix.index
	merged_tsne_df = TSNE_df.merge(node_df, on="strain")
	
	merged_tsne_df.to_csv('../Dataframes/TSNE_dataframe' + args.disease_name + '.csv', index = False)
	
	violin_plot_and_scatterplot_code(merged_tsne_df, "TSNE1", "TSNE2", "TSNE")
	
def UMAP():
	reducer = umap.UMAP(n_neighbors=200,
        min_dist=.05,
        n_components=2,
        init="spectral")
	X_transformed_umap = reducer.fit_transform(similarity_matrix)
	UMAP_df = pd.DataFrame(X_transformed_umap,columns=['UMAP' + str(i) for i in range(1,3)])
	UMAP_df["strain"] = similarity_matrix.index
	merged_umap_df = UMAP_df.merge(node_df, on="strain")
	
	merged_umap_df.to_csv('../Dataframes/UMAP_dataframe' + args.disease_name + '.csv', index = False)
	
	violin_plot_and_scatterplot_code(merged_umap_df, "UMAP1", "UMAP2", "UMAP")
	
def violin_plot_and_scatterplot_code(df, column1, column2, type_of_embedding)

	#creating violin plot
	violin_df = get_euclidean_data_frame(df, column1, column2, type_of_embedding)
	
	violin_df.to_csv('../Dataframes/' + type_of_embedding + '_violin_dataframe' + args.disease_name + '.csv', index = False)
	
	#creating scatterplot
	total_df = scatterplot_xyvalues(similarity_matrix.index, similarity_matrix, df, column1, column2, type_of_embedding)
	y_values = statsmodels.nonparametric.smoothers_lowess.lowess(
		total_df["euclidean"],
		total_df["genetic"],
		frac=0.6666666666666666,
		it=3,
		delta=0.0,
		is_sorted=False,
		missing='drop',
		return_sorted=True
	)

	PD_Y_values = pd.DataFrame(y_values)
	PD_Y_values.columns = ["LOWESS_x", "LOWESS_y"]

	total_df.to_csv('../Dataframes/' + type_of_embedding + '_xyvalues_scatterplot' + args.disease_name + '.csv', index = False)
	PD_Y_values.to_csv('../Dataframes/' + type_of_embedding + '_Y_values_scatterplot' + args.disease_name + '.csv', index = False)