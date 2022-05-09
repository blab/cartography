#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import linregress
import seaborn as sns


if __name__ == '__main__':
    # Configure command line interface.
    parser = argparse.ArgumentParser()
    parser.add_argument("--x", required=True, help="data frame of distances to plot on the x-axis")
    parser.add_argument("--y", required=True, help="data frame of distances to plot on the y-axis")
    parser.add_argument("--x-axis-label", required=True, help="x-axis label")
    parser.add_argument("--y-axis-label", required=True, help="y-axis label")
    parser.add_argument("--output", required=True, help="figure of distances plotted against each other")
    args = parser.parse_args()

    # Load distances.
    x = pd.read_csv(
        args.x,
        index_col=0,
    )

    y = pd.read_csv(
        args.y,
        index_col=0,
    )

    assert np.array_equal(x.index.values, y.index.values)

    x_array = squareform(x.values)
    y_array = squareform(y.values)

    regression = linregress(x_array, y_array)
    slope, intercept, r_value, p_value, std_err = regression

    fig, ax = plt.subplots(1, 1, figsize=(6, 6), dpi=200)
    ax.plot(
        x_array,
        y_array,
        "o",
        alpha=0.25,
    )

    ax.text(
        0.05,
        0.95,
        f"$R^2={r_value**2:.3f}$",
        horizontalalignment='left',
        verticalalignment='center',
        transform=ax.transAxes,
    )

    ax.set_xlabel(args.x_axis_label)
    ax.set_ylabel(args.y_axis_label)

    sns.despine()

    plt.tight_layout()
    plt.savefig(args.output)
