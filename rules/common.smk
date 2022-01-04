def _get_embedding_columns_by_wildcards(wildcards):
    method = wildcards.method.replace("-", "")

    if method in ("pca", "mds"):
        return f"{method}1 {method}2 {method}3 {method}4"
    else:
        return f"{method}_x {method}_y"
