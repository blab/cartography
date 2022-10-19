"""
takes the output from bases_missing.py and outputs a distance matrix
"""
import sys
sys.path.append("../")

import argparse
from augur.io import read_sequences
from collections import OrderedDict
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.spatial.distance import squareform, pdist

from Helpers import get_hamming_distances

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--alignment", required=True, help="aligned FASTA file of diseases")
    parser.add_argument("--output", required=True, help="path for the csv output")
    parser.add_argument("--indel-distance", action="store_true", help="add indel distances to genetic distances")

    args = parser.parse_args()

    sequences_by_name = OrderedDict()

    for sequence in read_sequences(args.alignment):
        sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names = list(sequences_by_name.keys())
    # Calculate Distance Matrix
    hamming_distances = get_hamming_distances(
        sequences_by_name.values(),
        args.indel_distance
    )
    distance_matrix = squareform(hamming_distances)

    distance_matrix = pd.DataFrame(distance_matrix, index=sequence_names)

    pd.DataFrame(distance_matrix).to_csv(args.output, index=True)
