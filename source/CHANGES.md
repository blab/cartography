# CHANGELOG

## 0.1.0

### Features

* Only calculate the distance matrix for an alignment if it isn't available already ([194fd74](https://github.com/blab/cartography/commit/194fd746c458d51bb73c962728da6c242a2d00f0))
* Source embedding params from cluster data ([8c26898](https://github.com/blab/cartography/commit/8c268981fa20d59888a92c0f38eedf8b42065db8))
* Initialize t-SNE embeddings with PCA instead of a random initialization ([89fc458](https://github.com/blab/cartography/commit/89fc4583e9af3caab332405dbb3b1f9e10f06c29))

### Bug Fixes

* Avoid re-reading alignment for PCA and t-SNE ([5fb2bbb](https://github.com/blab/cartography/commit/5fb2bbb13d686660dae48d1744f42a69c18fea57))

## 0.0.2

### Bug Fixes

* Issue [on github](https://github.com/blab/cartography/issues/20). The hamming distance calculation
now also works on lowercase fasta files as well as uppercase.


## 0.0.1

* First version of embed-pathogen.
