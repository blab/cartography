"""
Takes aligned strain data and outputs a fasta file that removes all low quality strains missing more than 3 standard deviations above the mean of missing bases"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from pathlib import Path
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, help="FASTA file of sequences")
    parser.add_argument("--output-fasta", required=True, help="name of disease")
    parser.add_argument("--output-bases-missing", help="the name of the output file for the PCA vs missing bases figure")
    parser.add_argument("--bases-missing-vs-pca", help="the path of the embed_pca file to output a missing bases plot of pca1 to missing bases")
        
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

    threshold = np.mean(genomes_missing_bases) + 3*np.std(genomes_missing_bases)

    new_strains = np.take(strains, np.where(genomes_missing_bases < threshold)).tolist()[0]
    new_genomes = np.take(genomes, np.where(genomes_missing_bases < threshold)).tolist()[0]

    if args.bases_missing_vs_pca is not None:

        genomes_missing_bases_df = pd.DataFrame(genomes_missing_bases, columns=["bases_missing"])
        genomes_missing_bases_df = genomes_missing_bases_df.merge(pd.DataFrame(strains, columns=["strain"]), how='outer', left_index = True, right_index = True)
        new_merged = genomes_missing_bases_df.merge(pd.read_csv(args.bases_missing_vs_pca, index_col=0), on="strain")
        
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        ax.plot(new_merged["bases_missing"], new_merged["pca1"], "o", alpha=0.25)

        ax.set_xlabel("Bases Missing")
        ax.set_ylabel(f"Euclidean distance (PCA1)")
        ax.set_title(f"Euclidean distance (PCA1) vs. Bases Missing")

        plt.savefig(args.output_bases_missing, dpi=300)

    strains_genomes = dict(zip(new_strains,new_genomes))

    if args.output_fasta is not None:
        fasta_file = open(args.output_fasta, "w")

        for i in range(len(new_strains)):

            fasta_file.write(">" + new_strains[i] + "\n" + new_genomes[i] + "\n")

        fasta_file.close()
        #with open(args.output_fasta, "w") as output_handle:
        #    SeqIO.write(strains_genomes, output_handle, "fasta")

