#!/usr/bin/env python3
import argparse
from augur.io import read_metadata
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="Nextstrain metadata with 'strain' and 'Nextclade_pango' fields")
    parser.add_argument("--priority-offset", type=float, default=0.25, help="offset for priorities of recombinant lineages relative to other lineages")
    parser.add_argument("--seed", type=int, default=314159, help="seed for random number generator")
    parser.add_argument("--output", required=True, help="tab-delimited floating point priorities per strain without a header")

    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    with open(args.output, "w", encoding="utf-8") as oh:
        metadata_reader = read_metadata(args.metadata, chunk_size=50000)
        for metadata in metadata_reader:
            metadata["priority"] = rng.random(size=metadata.shape[0])
            is_recombinant = metadata["Nextclade_pango"].str.startswith("X")
            metadata.loc[is_recombinant, "priority"] = metadata.loc[is_recombinant, "priority"] + args.priority_offset

            for strain, row in metadata.loc[:, ["priority"]].iterrows():
                print(f"{strain}\t{row['priority']:.2f}", file=oh)
