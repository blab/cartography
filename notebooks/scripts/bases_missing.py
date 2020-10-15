"""
Takes aligned strain data and outputs a two dataframes: strains_df and genomes_df, which can be converted to lists for hamming distance
"""
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform, pdist
from Bio import SeqIO
import re
from pathlib import Path
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, help="FASTA file of sequences")
    parser.add_argument("--output-fasta", required=True, help="name of disease")
        
    args = parser.parse_args()

    strains = []
    genomes = []
    for record in SeqIO.parse(args.alignment, "fasta"):
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
    genomes_missing_bases_df.columns = ["bases_missing", "strain"]

    threshold = np.mean(genomes_missing_bases) + 3*np.std(genomes_missing_bases)
    genomes_missing_bases_df = genomes_missing_bases_df[genomes_missing_bases_df['bases_missing'] >= threshold]
    list_of_val = genomes_missing_bases_df.index.tolist()

    strains = np.array(strains)
    genomes = np.array(genomes)

    deleted_strains = list(np.delete(strains, list_of_val))
    deleted_genomes = list(np.delete(genomes, list_of_val))

    if args.output_fasta is not None:
        fasta_file = open(args.output_fasta, "w")

        for i in range(len(deleted_genomes)):

            fasta_file.write(">" + deleted_strains[i] + "\n" + deleted_genomes[i] + "\n")

        fasta_file.close()

