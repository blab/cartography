"""Create time series cross-validation train/test splits and output data for a
given split.

"""
import argparse
from augur.io import read_sequences, write_sequences
import pandas as pd
from sklearn.model_selection import TimeSeriesSplit


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="metadata TSV to use for subsetting by time column")
    parser.add_argument("--alignment", required=True, help="alignment FASTA to subset")
    parser.add_argument("--distance-matrix", required=True, help="genetic distance matrix to subset")
    parser.add_argument("--time-column", required=True, help="integer column in the metadata to use for time series subsetting (e.g., generation, year, etc.)")
    parser.add_argument("--total-train-test-splits", type=int, required=True, help="total number of train/test splits to create")
    parser.add_argument("--train-test-split", type=int, required=True, help="index of specific train/test splits to output subset data for")
    parser.add_argument("--output-training-alignment", required=True, help="alignment FASTA for training")
    parser.add_argument("--output-training-genetic-distances", required=True, help="genetic distance matrix for training")
    parser.add_argument("--output-test-genetic-distances", required=True, help="genetic distance matrix for testing")

    args = parser.parse_args()

    # Load distances.
    distances = pd.read_csv(args.distance_matrix, index_col=0)
    distances.columns = distances.index.values

    # Load metadata.
    metadata = pd.read_csv(
        args.metadata,
        index_col=0,
        sep="\t",
    )

    # Filter metadata to only those strains considered in the distance matrix.
    metadata = metadata.loc[distances.index.values].copy()

    strains_by_time = {
        strain: set(strain_group.index.tolist())
        for strain, strain_group in metadata.groupby(args.time_column)
    }

    time_units = metadata[args.time_column].drop_duplicates().values
    time_series_cv = TimeSeriesSplit(n_splits=args.total_train_test_splits)

    for i, (train_time_index, test_time_index) in enumerate(time_series_cv.split(time_units)):
        # Only output data for the requested train/test split.
        if i == args.train_test_split:
            strains_train = set()
            for time_index in train_time_index:
                strains_train.update(strains_by_time[time_units[time_index]])

            strains_test = set()
            for time_index in test_time_index:
                strains_test.update(strains_by_time[time_units[time_index]])

            break

    # Only output training subset of the alignment.
    with open(args.output_training_alignment, "w") as oh:
        for sequence in read_sequences(args.alignment):
            if sequence.name in strains_train:
                write_sequences(sequence, oh)

    # Output training subset
    distances_train = distances.loc[strains_train, strains_train].copy()
    distances_train.to_csv(
        args.output_training_genetic_distances,
    )

    # Output test subset.
    distances_test = distances.loc[strains_test, strains_test].copy()
    distances_test.to_csv(
        args.output_test_genetic_distances,
    )
