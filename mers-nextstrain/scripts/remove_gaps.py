import argparse
import Bio.SeqIO
from collections import OrderedDict
import numpy as np
import pandas as pd
import re


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", help="a fasta file")
    parser.add_argument("--gaps", type=int, help="value of gaps per site to remove")
    parser.add_argument("--output-dataframe", help="a csv file")
    parser.add_argument("--output-fasta", help="a fasta file output")

    args = parser.parse_args()

    sequences_by_name = OrderedDict()

    for sequence in Bio.SeqIO.parse(args.alignment, "fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = list(sequences_by_name.keys())

    numbers = list(sequences_by_name.values())[:]
    for i in range(0,len(list(sequences_by_name.values()))):
        numbers[i] = re.sub(r'[^AGCT]', '5', numbers[i])
        numbers[i] = list(numbers[i].replace('A','1').replace('G','2').replace('C', '3').replace('T','4'))
        numbers[i] = [int(j) for j in numbers[i]]

    genomes_df = pd.DataFrame(numbers)
    genomes_df.columns = ["Site " + str(k) for k in range(0,len(numbers[i]))]

    # Drop columns with more than (number user picks) of bases missing in given site.
   
    dist_df = genomes_df.apply(pd.value_counts).iloc[4]

    dist_df = pd.DataFrame(dist_df)
    dist_df.columns = ["gaps"]

    print(dist_df)

    index_to_drop = dist_df.where(dist_df["gaps"] > args.gaps).dropna().index.values.tolist()
    print(index_to_drop)
    genomes_df = genomes_df.drop(index_to_drop, axis=1)

    if args.output_dataframe is not None:
        genomes_df.to_csv(args.output_dataframe)

    genomes = genomes_df.values.tolist()

    #converting to FASTA
    new_genomes = []
    if args.output_fasta is not None:
        for i in range(0, len(genomes)):
            string_numbers = ''.join(map(str, genomes[i])) 
            string_numbers = string_numbers.replace('1','A').replace('2','G').replace('3', 'C').replace('4','T').replace('5', "N")

            new_genomes.append(string_numbers)
    
        fasta_file = open(args.output_fasta, "w")

        for i in range(len(new_genomes)):

            fasta_file.write(">" + sequence_names[i] + "\n" + new_genomes[i] + "\n")

        fasta_file.close()