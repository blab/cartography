"""
Takes aligned strain data and PC1 and outputs a chart of bases missing vs PC1 
"""

import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, help="FASTA file of sequences")
    parser.add_argument("--pca", required=True, help="the path of the embed_pca file to output a missing bases plot of pca1 to missing bases")
    parser.add_argument("--output", help="the name of the output file for the PCA vs missing bases figure")
    
        
    args = parser.parse_args()

    strains = []
    genomes = []
    for record in SeqIO.parse(args.alignment, "fasta"):
        strains.append(str(record.id))
        genomes.append(str(record.seq))

    #Checking missing_bases

    genomes_missing_bases = []
    for x in genomes:
        x = re.sub(r'[^AGCT]', '5', x)
        numberOfN = x.count("5") #I'm leaving this logic here because MERS uses both "N" and "-" to dileneate missing sequences.
        genomes_missing_bases.append(numberOfN)


    genomes_missing_bases_df = pd.DataFrame(genomes_missing_bases, columns=["bases_missing"])
    genomes_missing_bases_df = genomes_missing_bases_df.merge(pd.DataFrame(strains, columns=["strain"]), how='outer', left_index = True, right_index = True)
    new_merged = genomes_missing_bases_df.merge(pd.read_csv(args.pca, index_col=0), on="strain")
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.plot(new_merged["bases_missing"], new_merged["pca1"], "o", alpha=0.25)

    ax.set_xlabel("Bases Missing")
    ax.set_ylabel(f"Euclidean distance (PCA1)")
    ax.set_title(f"Euclidean distance (PCA1) vs. Bases Missing")

    plt.savefig(args.output, dpi=300)

