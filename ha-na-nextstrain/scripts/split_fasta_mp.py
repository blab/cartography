import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", help="sequences to split into ha, na, and mp")
    parser.add_argument("--output_fastas", nargs=3, help="FASTA files of split genomes")

    args = parser.parse_args()

    strains = []
    genomes = []
    for record in SeqIO.parse(args.sequence, "fasta"):
        strains.append(str(record.id))
        genomes.append(str(record.seq))


    indices_na = []
    indices_ha = []
    indices_mp = []
    # Split genomes into two based on segment name
    for i in range(0, len(strains)):
        if(any(ele in strains[i] for ele in ["|N|4", "|S|4","|T|4"])):
            indices_ha.append(i)
        elif(any(ele in strains[i] for ele in ["|N|6", "|S|6","|T|6"])):
            indices_na.append(i)
        elif(any(ele in strains[i] for ele in ["|N|7", "|S|7","|T|7"])):
            indices_mp.append(i)
        else:
            print(strains[i])
            print("did not append anything")

    strains_na = np.array(strains)[indices_na]
    genomes_na = np.array(genomes)[indices_na]

    strains_mp = np.array(strains)[indices_mp]
    genomes_mp = np.array(genomes)[indices_mp]

    strains_ha = np.array(strains)[indices_ha]
    genomes_ha = np.array(genomes)[indices_ha]

    dictionary_ha = dict(zip(strains_ha, genomes_ha))
    dictionary_na = dict(zip(strains_na, genomes_na))
    dictionary_mp = dict(zip(strains_mp, genomes_mp))

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

    with open(args.output_fastas[2], "w") as output_handle:
        for strains, genomes in dictionary_mp.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")