import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", nargs="+", help="sequences to be concatenated together")
    parser.add_argument("--output_fasta", help="FASTA file of concatenated genomes")

    args = parser.parse_args()

    strains_ha = []
    genomes_ha = []
    for record in SeqIO.parse(args.sequences[0], "fasta"):
        strains_ha.append(str(record.id))
        genomes_ha.append(str(record.seq))
            
    strains_na = []
    genomes_na = []
    for record in SeqIO.parse(args.sequences[1], "fasta"):
        strains_na.append(str(record.id))
        genomes_na.append(str(record.seq))

    strings_matching = []
    strains = []
    for i in range(0,len(strains_ha)):
        if strains_ha[i] in strains_na:
            match = strains_na.index(strains_ha[i])
            strings_matching.append([i, match])
            strains.append(strains_ha[i])

    total_genome = []
    for i,j in strings_matching:
        concatenated_genome = genomes_ha[i] + genomes_na[j]
        total_genome.append(concatenated_genome)

    dictionary = dict(zip(strains, total_genome))

    #Output FASTA file
    with open(args.output_fasta, "w") as output_handle:
        for strains, genomes in dictionary.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")