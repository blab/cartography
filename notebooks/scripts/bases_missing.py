"""
Takes aligned strain data from snp sites and outputs a two dataframes: strains_df and genomes_df, which can be converted to lists for hamming distance
"""
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform, pdist
from Bio import SeqIO
import re
from pathlib import Path
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--sequences", required=True, help="FASTA file of sequences")
parser.add_argument("--dropped_strains", required=True, type=list, help="strains to drop")
parser.add_argument("--disease_name", required=True, help="name of disease")
	
args = parser.parse_args()

f = open(args.disease_name, "r")
virus_name = f.read()

dropped_strains = [line.rstrip('\n') for line in open(args.dropped_strains,"r")]

strains = []
genomes = []
for record in SeqIO.parse(sequences, "fasta"):
    if(record.id not in dropped_strains):
        strains.append(str(record.id))
        genomes.append(str(record.seq))
		
genomes_df = pd.DataFrame(genomes)
strains_df = pd.DataFrame(strains)
strains_df.columns = ["strain"]

#Checking missing_bases

genomes_missing_bases = []
for x in genomes:
    x = re.sub(r'[^AGCT]', '5', x)
    numberOfN = x.count("5")
    genomes_missing_bases.append(numberOfN)
	
genomes_missing_bases_df = pd.DataFrame(genomes_missing_bases)
genomes_missing_bases_df = genomes_missing_bases_df.merge(strains_df, how='outer', left_index = True, right_index = True)
genomes_missing_bases_df.columns = ["bases missing", "strain"]

#outputs a dataframe that can be visualized in the jupyter notebook - a histogram detailing amount of missing bases.

genomes_missing_bases_df.to_csv('notebooks/Dataframes/missing_bases' + virus_name + '.csv',index = False)

dropped_strains.extend(list(genomes_missing_bases_df[genomes_missing_bases_df["bases missing"]>1000]["strain"]))
dropped_strains = pd.Series(dropped_strains) 

# dropping all strains from the dataframes and lists that are in dropped_strains

indexNames = dropped_strains.isin(['strain']).index

strains_df.drop(indexNames , inplace=True)
strains_df.reset_index(drop = True)
strains_df.to_csv('notebooks/Dataframes/strains_dataframe' + virus_name + '.csv', index = False)

genomes_df.drop(indexNames , inplace=True)
genomes_df.columns = ["strain"]
genomes_df.reset_index(drop = True)
genomes_df.to_csv('notebooks/Dataframes/genomes_dataframe' + virus_name + '.csv', index = False)