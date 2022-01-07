"""
Takes aligned strain data and PC1 and outputs a chart of bases missing vs PC1
"""
import argparse
from augur.io import read_sequences
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from scipy.stats import linregress
import seaborn as sns

from Helpers import scatterplot_xyvalues
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, help="FASTA file of sequences")
    parser.add_argument("--embedding", required=True, help="the path of the embed_pca file to output a missing bases plot of pca1 to missing bases")
    parser.add_argument("--method", required=True, choices = ["pca", "mds", "t-sne", "umap"], help="the embedding used")
    parser.add_argument("--metadata", help="metadata if the color argument is true")
    parser.add_argument("--color", help="the column to color by (optional)")
    parser.add_argument("--bootstrapping-sample", default=100, type=int, help="number of times the data is sampled with replacement to find the mean and standard deviation of the pearson coefficient")
    parser.add_argument("--output", help="the name of the output file for the PCA vs missing bases figure")


    args = parser.parse_args()

    metadata_df = pd.read_csv(args.metadata, sep="\t")
    strains = []
    genomes = []
    for record in read_sequences(args.alignment):
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
    new_merged = genomes_missing_bases_df.merge(pd.read_csv(args.embedding, index_col=0), on="strain")

    r_value_arr = []
    for i in range(0, args.bootstrapping_sample):
        sampled_df = new_merged.sample(frac=1.0, replace=True)
        regression = linregress(new_merged["pca1"], new_merged["bases_missing"])
        slope, intercept, r_value, p_value, std_err = regression
        r_value_arr.append(r_value ** 2)

    r_value_arr = np.array(r_value_arr)

    mean = np.mean(r_value_arr)
    std = np.std(r_value_arr)
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

    ax.text(
                0.05,
                0.95,
                f"$R^2={mean:.3f} +/- {std}$",
                horizontalalignment='left',
                verticalalignment='center',
                transform=ax.transAxes,
            )

    plt.savefig(args.output, dpi=300)
