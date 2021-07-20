# Htool benchmarks

This repository contains some performance tests used by the library [Htool](https://htool-documentation.readthedocs.io/en/latest/). In particular, it is used to perform benchmarks with other libraries [here](https://github.com/PierreMarchand20/HMatrix_Benchmarks).

**Results can be find ...**

## General

### Geometries and kernels

The geometry and the kernels are defined [here](https://github.com/PierreMarchand20/HMatrix_Benchmarks#general-instructions). The following tests are performed for both Laplace and Helmholtz kernels, with both single- and double-precision. We use the following letter code to distinguish them:

- `S`: Laplace kernel, single precision,
- `D`: Laplace kernel, double precision,
- `C`: Helmholtz kernel, single precision,
- `Z`: Helmholtz kernel, double precision.

### Vectorization

- `none`: simple C++ code,
- `vcl`: vectorized with [Vector Class Library](https://github.com/vectorclass/version2),
- `xsimd`: vectorized with[xsimd](https://github.com/xtensor-stack/xsimd).

### Compressors

- `SVD`: Singular Value decomposition
- `fullACA`: full Adaptive Cross Approximation
- `partialACA`: partial Adaptive Cross Approximation
- `sympartialACA`: symmetric partial Adaptive Cross Approximation

> **Remark**
> In the future, we would like to add other compression techniques:
> 
> - Block ACA ([arxiv](https://arxiv.org/abs/1901.06101))
> - Generalized Nyström ([arxiv](https://arxiv.org/abs/2009.11392))


### Clustering techniques

- `PCARegularClustering`: axis computed via principal component analysis, split uniformly,
- `PCAGeometricClustering`: axis computed via principal component analysis, split barycentrically,
- `BoundingBox1RegularClustering`: axis along the largest extent of the bounding box, split uniformly,
- `BoundingBox1GeometricClustering`: axis along the largest extent of the bounding box, split barycentrically.

## Compression benchmark

We test all the previous compression techniques with all the different clustering techniques.

The output values are the rank of the low-rank approximation, assembly time, and matrix-vector product time.

## HMatrix benchmark

We only use `sympartialACA` with `PCARegularClustering`, but we use a different number of threads and MPI processes to test scalability.

The output values are the space-saving, assembly time, and the timing of 25 matrix-vector products.
