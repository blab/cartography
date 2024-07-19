#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import linregress


if __name__ == '__main__':
    # Configure command line interface.
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-genetic-distances", required=True, help="genetic distances for training data")
    parser.add_argument("--training-embedding", required=True, help="embedding based on training data only")
    parser.add_argument("--test-genetic-distances", required=True, help="genetic distances for test data")
    parser.add_argument("--test-embedding", required=True, help="embedding based on test data only")
    parser.add_argument("--output", help="CSV table with accuracy of embedding based on mean absolute error of observed and estimated test data")
    parser.add_argument("--annotations", nargs="*", help="annotations to add to statistics table")
    args = parser.parse_args()

    # Load distances.
    training_genetic_distances = pd.read_csv(
        args.training_genetic_distances,
        index_col=0,
    ).sort_index(axis=0).sort_index(axis=1)

    test_genetic_distances = pd.read_csv(
        args.test_genetic_distances,
        index_col=0,
    ).sort_index(axis=0).sort_index(axis=1)

    # Load training embedding.
    training_embedding = pd.read_csv(
        args.training_embedding,
        index_col=0,
    ).sort_index(axis=0)

    # Load test embedding.
    test_embedding = pd.read_csv(
        args.test_embedding,
        index_col=0,
    ).sort_index(axis=0)

    # Fit a linear model between genetic and embedding distances.
    training_embedding_distances = pdist(training_embedding.values)
    training_genetic_distances = squareform(training_genetic_distances)
    model = linregress(training_embedding_distances, training_genetic_distances)

    # Estimate genetic distances from test embedding distances.
    test_embedding_distances = pdist(test_embedding.values)
    observed_test_genetic_distances = squareform(test_genetic_distances)
    estimated_test_genetic_distances = model.intercept + (model.slope * test_embedding_distances)

    # Calculate mean absolute error (MAE) between observed and estimated
    # genetic distances for test data.
    error = (np.abs(observed_test_genetic_distances - estimated_test_genetic_distances)).sum() / observed_test_genetic_distances.shape[0]

    statistics = pd.DataFrame([{
        "mae": error,
    }])

    if args.annotations:
        for annotations in args.annotations:
            annotations = eval(annotations)
            for key, value in annotations.items():
                statistics[key] = value

    statistics.to_csv(
        args.output,
        index=False,
    )
