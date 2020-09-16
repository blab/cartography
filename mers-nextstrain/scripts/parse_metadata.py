"""Parse metadata from a FASTA headers like augur parse but only output the metadata table and use the FASTA header as the strain name.
"""
import argparse
from augur.utils import numeric_date
import Bio.SeqIO
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse metadata from a FASTA headers like augur parse but only output the metadata table and use the FASTA header as the strain name.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--sequences", help="a FASTA of sequences with |-delimited metadata in the defline")
    parser.add_argument("--fields", nargs="+", help="names for fields in the given FASTA header")
    parser.add_argument("--delimiter", default="|", help="delimiter for fields in the FASTA header")
    parser.add_argument("--output-sequences", help="a FASTA field of sequences")
    parser.add_argument("--output-metadata", help="a metadata TSV file")

    args = parser.parse_args()

    sequences = Bio.SeqIO.parse(args.sequences, "fasta")

    cleaned_sequences = []
    metadata_records = []
    for sequence in sequences:
        # Jeddah-1|KF917527|camel|2013-11-08 -> Jeddah-1-KF917527-camel-2013-11-08
        sequence_id = sequence.id.replace(args.delimiter, "-")
        attributes = sequence.id.split(args.delimiter)
        metadata_record = dict(zip(args.fields, attributes))

        if "strain" in metadata_record:
            metadata_record["original_strain"] = metadata_record["strain"]

        metadata_record["strain"] = sequence_id
        metadata_records.append(metadata_record)

        sequence.id = sequence_id
        sequence.description = ""
        cleaned_sequences.append(sequence)

    df = pd.DataFrame(metadata_records)

    # Annotate a numerical date.
    if "date" in args.fields:
        df["num_date"] = pd.to_datetime(df["date"]).apply(numeric_date)

    # Save the metadata.
    df.to_csv(args.output_metadata, sep="\t", index=False)

    # Save the sequences.
    Bio.SeqIO.write(cleaned_sequences, args.output_sequences, "fasta")
