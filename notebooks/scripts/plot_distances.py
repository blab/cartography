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
    parser.add_argument("--output-figure", required=True, help="figure of distances plotted against each other")
    parser.add_argument("--output-statistics", help="table of statistics from the linear regression for the given plot")
    parser.add_argument("--annotations", nargs="*", help="annotations to add to statistics table")
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

    assert np.array_equal(sorted(x.index.values), sorted(y.index.values))
    x_indices_in_y = [np.where(y.index.values == x_value)[0][0] for x_value in x.index.values]
    x_values = x.values
    y_values = y.iloc[x_indices_in_y, x_indices_in_y].values

    x_array = squareform(x_values)
    y_array = squareform(y_values)

    regression = linregress(x_array, y_array)
    slope, intercept, r_value, p_value, std_err = regression

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=200)
    # ax.plot(
    #     x_array,
    #     y_array,
    #     "o",
    #     alpha=0.25,
    # )
    data = pd.DataFrame({"x": x_array, "y": y_array})
    # sns.boxplot(
    #     x="x",
    #     y="y",
    #     data=data,
    #     linewidth=0.75,
    #     fliersize=0.1,
    #     color="#CCCCCC",
    #     ax=ax
    # )
    g = sns.jointplot(
        x=x_array,
        y=y_array,
        kind="hex",
        color="#4CB391",
        height=8,
    )
    g.set_axis_labels(args.x_axis_label, args.y_axis_label)

    g.ax_joint.text(
        0.05,
        0.95,
        f"$R^2={r_value**2:.3f}$",
        horizontalalignment='left',
        verticalalignment='center',
        transform=g.ax_joint.transAxes,
    )

    #ax.set_xlabel(args.x_axis_label)
    #ax.set_ylabel(args.y_axis_label)

    sns.despine()

    #plt.tight_layout()
    plt.savefig(args.output_figure)

    if args.output_statistics:
        statistics = pd.DataFrame([{
            "slope": slope,
            "intercept": intercept,
            "r_value": r_value,
            "r_squared": r_value ** 2,
            "p_value": p_value,
            "std_err": std_err,
        }])

        if args.annotations:
            for annotations in args.annotations:
                annotations = eval(annotations)
                for key, value in annotations.items():
                    statistics[key] = value

        statistics.to_csv(
            args.output_statistics,
            index=False,
        )
