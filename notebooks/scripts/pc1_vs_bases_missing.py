"""
Takes aligned strain data and PC1 and outputs a chart of bases missing vs PC1 
"""

import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, help="FASTA file of sequences")
    parser.add_argument("--metadata", help="metadata if the color argument is true")
    parser.add_argument("--pca", required=True, help="the path of the embed_pca file to output a missing bases plot of pca1 to missing bases")
    parser.add_argument("--color", help="the column to color by (optional)")
    parser.add_argument("--output", help="the name of the output file for the PCA vs missing bases figure")
    
        
    args = parser.parse_args()

    metadata_df = pd.read_csv(args.metadata, sep="\t")
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
    if args.color is not None:
        new_merged = new_merged.merge(metadata_df, on="strain")
        
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    if args.color is None:
        sns.scatterplot(data=new_merged, x="bases_missing", y="pca1", alpha=0.25, ax=ax)
    else:
        sns.scatterplot(data=new_merged, x="bases_missing", y="pca1", hue=args.color, alpha=0.25)
    if args.color is not None:
        ax.legend(title="low quality", loc="upper left")

    ax.set_xlabel("Bases Missing")
    ax.set_ylabel(f"PC1")

    plt.savefig(args.output, dpi=300)

