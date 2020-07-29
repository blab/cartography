"""
takes the output from bases_missing.py and outputs a distance matrix
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.spatial.distance import squareform, pdist
from Helpers import get_hamming_distances


f = open(snakemake.input.disease_name, "r")
virus_name = f.read()

strains_df = pd.read_csv('notebooks/Dataframes/strains_dataframe.csv')  
genomes_df = pd.read_csv('notebooks/Dataframes/genomes_dataframe.csv')

genomes = genomes_df["strain"].values.tolist()

strains = strains_df["strain"].values.tolist()

%%time

# Calculate Hamming distances.
hamming_distances = get_hamming_distances(genomes)

# Convert distinct pairwise distances into the more redundant but more interpretable square matrix.
similarity_matrix = squareform(hamming_distances)

# Convert the numpy matrix to a pandas data frame with strain annotations for rows and columns.
similarity_matrix = pd.DataFrame(
    similarity_matrix,
    columns=strains,
    index=strains
)

# Write out the resulting data frame to cache distance calculations.
# Keep the index in the output file, so it is immediately available on read.
similarity_matrix.to_csv('notebooks/Dataframes/distance_matrix_' + virus_name + '.csv')