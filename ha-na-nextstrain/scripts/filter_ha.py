import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", nargs=2, help="sequences to filter, the one to filter first and one to filter by second")
    parser.add_argument("--output_fasta", nargs=2, help="FASTA files of split genomes")

    args = parser.parse_args()

    strains_filter = []
    genomes_filter = []
    for record in SeqIO.parse(args.sequence[0], "fasta"):
        strains_filter.append(str(record.id))
        genomes_filter.append(str(record.seq))

    strains_filter_by = []
    genomes_filter_by = []
    for record in SeqIO.parse(args.sequence[1], "fasta"):
        strains_filter_by.append(str(record.id))
        genomes_filter_by.append(str(record.seq))

    strains_ = []
    for i in range(0, len(strains_filter)):
        if strains_filter[i] in strains_filter_by:
            strains_.append(i)

    strains_ha = np.array(strains_filter)[strains_]
    genomes_ha = np.array(genomes_filter)[strains_]

    dictionary_ha = dict(zip(strains_ha, genomes_ha))
    dictionary_na = dict(zip(strains_filter_by, genomes_filter_by))

    with open(args.output_fasta[0], "w") as output_handle:
        for strains, genomes in dictionary_ha.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")

    with open(args.output_fasta[1], "w") as output_handle:
        for strains, genomes in dictionary_na.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")