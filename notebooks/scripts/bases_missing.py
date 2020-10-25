"""
Takes aligned strain data and outputs a fasta file that removes all low quality strains missing more than 3 standard deviations above the mean of missing bases"""

import numpy as np
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

    #Checking missing_bases

    genomes_missing_bases = []
    for x in genomes:
        x = re.sub(r'[^AGCT]', '5', x)
        numberOfN = x.count("5") #I'm leaving this logic here because MERS uses both "N" and "-" to dileneate missing sequences.
        genomes_missing_bases.append(numberOfN)

    threshold = np.mean(genomes_missing_bases) + 3*np.std(genomes_missing_bases)

    new_strains = np.take(strains, np.where(genomes_missing_bases >= threshold)).tolist()[0]
    new_genomes = np.take(genomes, np.where(genomes_missing_bases >= threshold)).tolist()[0]

    strains_genomes = dict(zip(new_strains,new_genomes))

    if args.output_fasta is not None:
        with open(args.output_fasta, "w") as output_handle:
            SeqIO.write(strains_genomes, output_handle, "fasta")

