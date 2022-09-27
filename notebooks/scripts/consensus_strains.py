import argparse
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
from augur.io import read_sequences
from collections import OrderedDict
import pandas as pd
from augur.utils import read_node_data
import pathlib

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="the name of the aligned file to create consensus strains from")
    parser.add_argument("--metadata", required=True, help="clade information to create cluster per clade.")
    parser.add_argument("--group-by", default="clade_membership", help="attribute name in metadata to group sequences from the alignment to find a consensus for")
    parser.add_argument("--output", required=True, help="outputting a csv file of metadata info for HDBSCAN results")

    args = parser.parse_args()

    sequences_by_name = OrderedDict()

    for sequence in read_sequences(args.alignment):#"../ha-na-nextstrain/results/aligned_concatenated.fasta"):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = list(sequences_by_name.keys())

    if ("json" in pathlib.Path(args.metadata).suffix):

        node_data = read_node_data(args.metadata)

        clade_annotations = pd.DataFrame([
        {"strain": strain, "clade": annotations["MCC"]}
            for strain, annotations in node_data["nodes"].items()
            if strain in sequences_by_name
        ])
        clade_annotations = clade_annotations.set_index("strain")
    else:
        clade_annotations = pd.read_csv(args.metadata, sep="\t").set_index("strain")
        clade_annotations = clade_annotations[["clade_membership"]]
        clade_annotations.columns=["clade"]

    clade = clade_annotations.groupby(["clade"])
    list_clades = [clade.get_group(x) for x in clade.groups]
    # merge strain information (make a dictionary)
    for clades in list_clades:
        sequences = []
        records = []
        for index, row in clades.iterrows():
            sequences.append(sequences_by_name[index])
            records.append(
                SeqRecord(Seq(sequences_by_name[index]), id=index)
            )
        clades["sequence"] = sequences

        alignment = MultipleSeqAlignment(records)
        consensus_sequence = SummaryInfo(alignment).dumb_consensus(threshold=0.5)
        alignment.append(SeqRecord(consensus_sequence, id="consensus"))

        cluster_id = clades["clade"].drop_duplicates().values[0]
        with open(f"alignment_{cluster_id}.fasta", "w") as handle:
            SeqIO.write(alignment, handle, "fasta")

    total_consensus_per_group = []
    for clades in list_clades:
        consensus = []
        for (name, data) in clades[["sequence"]].iteritems():
            split_data = pd.Series(data[0]).apply(lambda x: pd.Series(list(x)))
            consensus.append(split_data.value_counts().idxmax())
        str_consensus = "".join(list(consensus[0]))
        total_consensus_per_group.append(str_consensus)

    clade_consensus_strain = dict(zip(clade.groups, total_consensus_per_group))

    # fasta_output = []
    # fasta_output.append()
    # for clade, consensus in zip(clade.groups, total_consensus_per_group):
    #     record = SeqRecord(consensus, clade)
    #     fasta_output.append(record)

    # SeqIO.write(fasta_output, args.output, "fasta")

    #Output FASTA file
    with open(args.output, "w") as output_handle:
        for strains, genomes in clade_consensus_strain.items():
            record = SeqRecord(
            Seq(genomes),
            id=strains,
            description=""
            )
            SeqIO.write(record, output_handle, "fasta")
