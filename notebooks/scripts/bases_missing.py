"""
Takes aligned strain data and outputs a fasta file that removes all low quality strains missing more than 3 standard deviations above the mean of missing bases"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
        numberOfN = x.count("5") #This logic is here because MERS uses both "N" and "-" to dileneate missing sequences.
        genomes_missing_bases.append(numberOfN)

    threshold = np.mean(genomes_missing_bases) + 3*np.std(genomes_missing_bases)

    #create dict
    dictionary = dict(zip(np.take(strains, np.where(genomes_missing_bases < threshold)).tolist()[0], np.take(genomes, np.where(genomes_missing_bases < threshold)).tolist()[0]))

    with open(args.output_fasta, "w") as output_handle:
        for strains, genomes in dictionary.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")