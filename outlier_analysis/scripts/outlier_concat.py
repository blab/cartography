import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import re


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-human", help="sequence - main FASTA")
    parser.add_argument("--sequence-outlier", help= "outlier sequences")
    parser.add_argument("--metadata-outlier", help= "outlier metadata")
    parser.add_argument("--metadata-human", help="metadata - main TSV")
    parser.add_argument("--outlier-list", required=False, help="list of outliers within the main file, if they exist")
    parser.add_argument("--output-fasta", help="FASTA file of concatenated genomes")
    parser.add_argument("--output-metadata", help="TSV file of metadata for genomes")

    args = parser.parse_args()


    # read in outliers
    with open(args.outlier_list, 'r') as outlier:
        lines = outlier.read().split('\n')

    # concatenate metadata
    human_df = pd.read_csv(args.metadata_human, sep='\t')
    outlier_df = pd.read_csv(args.metadata_outlier, sep='\t')
    concatenated_metadata = pd.concat([human_df, outlier_df])

    lines = lines + outlier_df["strain"].values.tolist()

    strain_names = concatenated_metadata['strain'].values.tolist()
    outlier_list = []

    for i in range(0,len(lines)):
        lines[i] = re.sub(r'\|.*$', '', str(lines[i]))  
        
    for i in range(0,len(strain_names)):
        if strain_names[i] in lines:
            outlier_list.append("outlier")
        else:
            outlier_list.append("not_outlier")

    concatenated_metadata['outlier'] = outlier_list


    # read out metadata
    concatenated_metadata.to_csv(args.output_metadata, index=False, sep="\t")

    # read in and concatenate sequences
    strains_human = []
    genomes_human = []
    for record in SeqIO.parse(args.sequence_human, "fasta"):
        strains_human.append(str(record.id))
        genomes_human.append(str(record.seq))
            
    strains_outlier = []
    genomes_outlier = []
    for record in SeqIO.parse(args.sequence_outlier, "fasta"):
        strains_outlier.append(str(record.id))
        genomes_outlier.append(str(record.seq))

    total_strains = strains_human + strains_outlier
    total_genomes = genomes_human + genomes_outlier
 
    dictionary = dict(zip(total_strains, total_genomes))

    #Output FASTA file
    with open(args.output_fasta, "w") as output_handle:
        for strains, genomes in dictionary.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")