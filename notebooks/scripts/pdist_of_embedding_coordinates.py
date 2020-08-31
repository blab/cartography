import argparse
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--embedding", required=True, help="aligned FASTA file of diseases")
    parser.add_argument("--output", required=True, help="path for the csv output")

    args = parser.parse_args()

    embedding_df = pd.read_csv(args.embedding, index_col=0)

    embedding_df["strain"] = embedding_df.index

    scaler = StandardScaler()

    scaled_distance = scaler.fit_transform(pdist(embedding_df.drop(["strain"], axis = 1)).reshape(-1, 1))

    pd.DataFrame(scaled_distance).to_csv(args.output)