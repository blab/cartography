import Bio.SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", help="gb file to convert")
    parser.add_argument("--output", help="fasta file output")
    args = parser.parse_args()


    record = next(Bio.SeqIO.parse(args.file, "genbank"))
    Bio.SeqIO.write(record, args.output, "fasta")







