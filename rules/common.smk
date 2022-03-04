def _get_embedding_columns_by_wildcards(wildcards):
    method = wildcards.method.replace("-", "")

    if method == "pca":
        return f"{method}1 {method}2 {method}3 {method}4 {method}5 {method}6 {method}7 {method}8 {method}9 {method}10"
    elif method == "mds":
        return f"{method}1 {method}2 {method}3 {method}4"
    else:
        return f"{method}_x {method}_y"
