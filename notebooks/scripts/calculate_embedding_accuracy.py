#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform
from scipy.stats import linregress


if __name__ == '__main__':
    # Configure command line interface.
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-embedding", required=True, help="full embedding with train and test strains included")
    parser.add_argument("--training-genetic-distances", required=True, help="genetic distances for training data")
    parser.add_argument("--test-genetic-distances", required=True, help="genetic distances for test data")
    parser.add_argument("--training-embedding", required=True, help="embedding based on training data only")
    parser.add_argument("--output", help="CSV table with accuracy of embedding based on mean squared error of observed and estimated test data")
    parser.add_argument("--annotations", nargs="*", help="annotations to add to statistics table")
    args = parser.parse_args()

    # Load full embedding.
    full_embedding = pd.read_csv(
        args.full_embedding,
        index_col=0,
    )

    # Load distances.
    training_genetic_distances = pd.read_csv(
        args.training_genetic_distances,
        index_col=0,
    )

    test_genetic_distances = pd.read_csv(
        args.test_genetic_distances,
        index_col=0,
    )

    # Load training embedding.
    training_embedding = pd.read_csv(
        args.training_embedding,
        index_col=0,
    )

    # Pad training embedding with zeros to have the same number of records as
    # the full embedding for Procrustes projection.
    padding_length = full_embedding.shape[0] - training_embedding.shape[0]
    padding_width = training_embedding.shape[1]
    padded_training_embedding = np.concatenate((training_embedding, np.zeros((padding_length, padding_width))))

    # Project training embedding into full embedding space.
    centered_full_embedding, padded_projected_training_embedding, disparity = procrustes(
        full_embedding,
        padded_training_embedding,
    )

    # Calculate pairwise distances from full embedding after standardization
    # with Procrustes.
    full_embedding_distances = squareform(pdist(centered_full_embedding))
    full_embedding_distances = pd.DataFrame(
        full_embedding_distances,
        index=full_embedding.index.values,
        columns=full_embedding.index.values,
    )

    # Remove padding from projected embedding.
    projected_training_embedding = padded_projected_training_embedding[:training_embedding.shape[0]]

    # Fit a linear model between genetic and embedding distances.
    training_genetic_distances = squareform(training_genetic_distances)
    training_embedding_distances = pdist(projected_training_embedding)
    model = linregress(training_genetic_distances, training_embedding_distances)

    # Estimate embedding distances from test genetic distances.
    strains_test = test_genetic_distances.index.tolist()
    test_genetic_distances = squareform(test_genetic_distances)
    test_embedding_distances = model.intercept + (model.slope * test_genetic_distances)

    # Calculate mean squared error (MSE) between observed and estimated
    # embedding distances for test data.
    known_test_embedding_distances = squareform(
        full_embedding_distances.loc[
            strains_test,
            strains_test
        ]
    )
    mse = ((known_test_embedding_distances - test_embedding_distances) ** 2).sum() / known_test_embedding_distances.shape[0]

    statistics = pd.DataFrame([{
        "mse": mse,
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
