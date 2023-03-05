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
    parser.add_argument("--time-column", required=True, help="integer column in the metadata to use for time series subsetting (e.g., generation, year, etc.)")
    parser.add_argument("--total-train-test-splits", type=int, required=True, help="total number of train/test splits to create")
    parser.add_argument("--train-test-split", type=int, required=True, help="index of specific train/test splits to output subset data for")
    parser.add_argument("--train-test-size", type=int, required=True, help="number of indices to include in train/test splits")
    parser.add_argument("--gap", default=0, type=int, help="gap between train/test time points in units of the given time column")
    parser.add_argument("--output-training-alignment", required=True, help="alignment FASTA for training")
    parser.add_argument("--output-test-alignment", required=True, help="alignment FASTA for testing")

    args = parser.parse_args()

    # Load metadata.
    metadata = pd.read_csv(
        args.metadata,
        index_col=0,
        sep="\t",
    )

    strains_by_time = {
        strain: set(strain_group.index.tolist())
        for strain, strain_group in metadata.groupby(args.time_column)
    }

    # Find all available time units (e.g., generations) and calculate a range of
    # continuous values between the minimum and maximum values. These values
    # will be the input to the time series split logic below.
    time_units = metadata[args.time_column].drop_duplicates().values
    continuous_time_units = list(range(time_units.min(), time_units.max() + 1))

    # Generate time series train/test splits with an equal size for both train
    # and test data.
    time_series_cv = TimeSeriesSplit(
        n_splits=args.total_train_test_splits,
        max_train_size=args.train_test_size,
        test_size=args.train_test_size,
        gap=args.gap,
    )

    for i, (train_time_index, test_time_index) in enumerate(time_series_cv.split(continuous_time_units)):
        # Only output data for the requested train/test split.
        if i == args.train_test_split:
            strains_train = set()
            for time_index in train_time_index:
                # Not all values in the continous range of time units will be
                # sampled by the simulation output, so check for whether the
                # current time index has samples or not.
                time_unit = continuous_time_units[time_index]
                if time_unit in strains_by_time:
                    strains_train.update(strains_by_time[time_unit])

            strains_test = set()
            for time_index in test_time_index:
                time_unit = continuous_time_units[time_index]
                if time_unit in strains_by_time:
                    strains_test.update(strains_by_time[time_unit])

            break

    # Only output training or test subsets of the alignment.
    with open(args.output_training_alignment, "w") as oh_training:
        with open(args.output_test_alignment, "w") as oh_test:
            for sequence in read_sequences(args.alignment):
                if sequence.name in strains_train:
                    write_sequences(sequence, oh_training)
                elif sequence.name in strains_test:
                    write_sequences(sequence, oh_test)
