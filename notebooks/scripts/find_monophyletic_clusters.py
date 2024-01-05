#!/usr/bin/env python3
import argparse
from augur.utils import read_node_data, read_tree
from collections import Counter
import csv


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dataset-name", required=True, help="name of the dataset being analyzed")
    parser.add_argument("--tree", required=True, help="Newick tree with named internal nodes (from augur refine)")
    parser.add_argument("--labels", required=True, help="node data JSON with inferred cluster labels for internal nodes in the given tree")
    parser.add_argument("--methods", required=True, nargs="+", help="list of methods that should match the list of columns with cluster labels to analyze")
    parser.add_argument("--attributes", required=True, nargs="+", help="list of node data attributes with cluster labels to analyze")
    parser.add_argument("--output", required=True, help="CSV output summarizing the number of clusters and cluster transitions per cluster label column")

    args = parser.parse_args()

    # Load tree.
    tree = read_tree(args.tree)

    # Load node data.
    node_data = read_node_data(args.labels)["nodes"]

    # Count the number of transitions between internal node states on the tree
    # for each method/attribute.
    with open(args.output, "w", encoding="utf-8") as oh:
        csv_writer = csv.writer(oh, lineterminator="\n")
        csv_writer.writerow((
            "dataset",
            "method",
            "n_clusters",
            "n_cluster_transitions",
            "n_extra_transitions",
            "clusters",
            "transitions",
        ))

        for method, label_attribute in zip(args.methods, args.attributes):
            distinct_clusters = set()
            cluster_transitions = []
            for node in tree.find_clades(order="preorder", terminal=False):
                node_cluster = node_data[node.name][label_attribute]
                distinct_clusters.add(node_cluster)

                for child in node.clades:
                    child_cluster = node_data[child.name][label_attribute]
                    distinct_clusters.add(child_cluster)

                    if not child.is_terminal() and child_cluster != node_cluster and child_cluster != "-1":
                        cluster_transitions.append((node_cluster, child_cluster))

            distinct_clusters = sorted(distinct_clusters)
            cluster_transition_counts = Counter(cluster_transitions)
            extra_transitions = sum(value - 1 for value in cluster_transition_counts.values() if value > 1)

            csv_writer.writerow((
                args.dataset_name,
                method,
                len(distinct_clusters),
                len(cluster_transitions),
                extra_transitions,
                f"{distinct_clusters!r}",
                f"{cluster_transitions!r}",
            ))
