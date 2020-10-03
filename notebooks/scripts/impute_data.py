# converts unimputed data into a genome_df, creates imputed data
import argparse
import Bio.SeqIO
from collections import OrderedDict
import numpy as np
import pandas as pd
import re
from sklearn.impute import SimpleImputer


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "creates imputed fasta or df", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        
    parser.add_argument("--alignment", help="the path to an alignment FASTA file")
    parser.add_argument("--output-dataframe", help="path for outputting as a dataframe of numbers")
    parser.add_argument("--output-fasta", help="path for outputting as a FASTA")

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
    #performing PCA on my pandas dataframe

    # Impute missing values using the most frequent value in each column.
    imputer = SimpleImputer(missing_values=5, strategy="most_frequent")
    genomes_df = pd.DataFrame(imputer.fit_transform(genomes_df))
    genomes_df.index = sequence_names
    if args.output_dataframe is not None:
        genomes_df.to_csv(args.output_dataframe)

    genomes = genomes_df.values.tolist()
    #converting to FASTA
    new_genomes = []
    if args.output_fasta is not None:
        for i in range(0, len(genomes)):
            string_numbers = ''.join(map(str, genomes[i])) 
            string_numbers = string_numbers.replace('1','A').replace('2','G').replace('3', 'C').replace('4','T')

            new_genomes.append(string_numbers)
    
        fasta_file = open(args.output_fasta, "w")

        for i in range(len(new_genomes)):

            fasta_file.write(">" + sequence_names[i] + "\n" + new_genomes[i] + "\n")

        fasta_file.close()