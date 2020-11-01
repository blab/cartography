import argparse
from dendropy.calculate import popgenstat
from dendropy import DnaCharacterMatrix

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "find nucleotide diversity of a population", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--alignment", help="an aligned FASTA file to create a DNA character matrix from")
    parser.add_argument("--output", help="outputting a txt file with the nucleotide_diversity value")

    args = parser.parse_args()

    d = DnaCharacterMatrix.get(
        path=args.alignment,
        schema="fasta"
    )

    with open(args.output, 'w') as f:
        f.write(str(popgenstat.nucleotide_diversity(d, ignore_uncertain=True)))