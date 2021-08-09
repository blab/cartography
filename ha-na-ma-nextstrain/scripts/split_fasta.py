import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", help="sequences to split into ha and na")
    parser.add_argument("--output_fastas", nargs=2, help="FASTA files of split genomes")

    args = parser.parse_args()

    strains = []
    genomes = []
    for record in SeqIO.parse(args.sequence, "fasta"):
        strains.append(str(record.id))
        genomes.append(str(record.seq))


    indices_na = []
    indices_ha = []
    # Split genomes into two based on segment name
    for i in range(0, len(strains)):
        if(any(ele in strains[i] for ele in ["|N|4", "|S|4","|T|4"])):
            indices_ha.append(i)
        else:
            indices_na.append(i)

    strains_na = np.array(strains)[indices_na]
    genomes_na = np.array(genomes)[indices_na]

    strains_ha = np.array(strains)[indices_ha]
    genomes_ha = np.array(genomes)[indices_ha]

    dictionary_ha = dict(zip(strains_ha, genomes_ha))
    dictionary_na = dict(zip(strains_na, genomes_na))

    #Output FASTA files
    with open(args.output_fastas[0], "w") as output_handle:
        for strains, genomes in dictionary_ha.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")

    with open(args.output_fastas[1], "w") as output_handle:
        for strains, genomes in dictionary_na.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")