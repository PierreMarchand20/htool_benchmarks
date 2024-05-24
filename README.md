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
- `xsimd`: vectorized with [xsimd](https://github.com/xtensor-stack/xsimd).

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

## Task-based implementation benchmark
### Overview
bench_hmatrix_build.cpp uses the Google Benchmark library to benchmark two different implementations of a function called `test_hmatrix_build`.

The program defines two benchmark functions: `BM_test_hmatrix_build` which use the classic implementation of the `HMatrixTreeBuilder` class and `BM_NEW_test_hmatrix_build` which use the task-based implementation of the `HMatrixTaskBasedTreeBuilder` class. These functions are called by the Google Benchmark library to measure the performance of the `test_hmatrix_build` function.

The benchmark functions take a benchmark::State object as a parameter, which provides information about the benchmark execution, such as the number of iterations and the current iteration number. Inside the benchmark functions, the `test_hmatrix_build` function is called multiple times, and the time and memory usage are measured using the state object.

The program also includes some constants and defines the benchmarks using the BENCHMARK macro provided by the Google Benchmark library. The benchmarks specify the range of input values, the number of repetitions, the minimum warm-up time, the number of threads, and the minimum time per repetition.

It is recommended to run bench_hmatrix_build.cpp with the command "./scripts/run_bench_treebuilder.sh" because this script prepares the machine for benchmarking by changing the CPU governor to performance mode, disabling turbo boost, and disabling ASLR (Address Space Layout Randomization). These commands aim at reducing the noise in the time and memory usage measurement. It then runs the benchmark program bench_hmatrix_build with specific options, such as writing the benchmark results to a file in JSON format, enabling random interleaving of benchmark repetitions, and setting the task affinity to a specific CPU. Then the results are compared by a python script "compare.py". After running the benchmark, it restores the machine settings back to powersave mode, enables turbo boost, and enables ASLR. 

### Interpreting the output 
A breakdown of each column:
    - `Benchmark`: The name of the function being benchmarked, along with the size of the input (after the slash).
    - `Time`: The average time per operation, across all iterations.
    - `CPU`: The average CPU time per operation, across all iterations.
    - `Iterations`: The number of iterations the benchmark was run to get a stable estimate.
    - `Time Old` and `Time New`: These represent the average time it takes for a function to run in two different scenarios or versions. For example, you might be comparing how fast a function runs before and after you make some changes to it.
    - `CPU Old` and `CPU New`: These show the average amount of CPU time that the function uses in two different scenarios or versions. This is similar to `Time Old` and `Time New`, but focuses on CPU usage instead of overall time.
In the comparing section the values in `Time` and `CPU` columns are calculated as (new - old) / |old|.

A breakdown of each row:
    - `*_pvalue`: This shows the p-value for the statistical test comparing the performance of the process running with one thread. A value of 0.0000 suggests a statistically significant difference in performance. The comparison was conducted using the U Test (Mann-Whitney U Test) with at least 9 repetitions for each case. The result of said the statistical test is additionally communicated through color coding. Green means that the benchmarks are statistically different, this could mean the performance has either significantly improved or significantly deteriorated. You should look at the actual performance numbers to see which is the case. Red means that the benchmarks are statistically similar, this means the performance hasn't significantly changed.
    - `*_mean`: This shows the relative difference in mean execution time between two different cases. 
    - `*_median`: Similarly, this shows the relative difference in the median execution time. 
    - `*_stddev`: This is the relative difference in the standard deviation of the execution time, which is a measure of how much variation or dispersion there is from the mean. 
    - `*_cv`: CV stands for Coefficient of Variation. It is the ratio of the standard deviation to the mean. It provides a standardized measure of dispersion. 
    - `OVERALL_GEOMEAN`: Geomean stands for geometric mean, a type of average that is less influenced by outliers. If the values are all zero for the old and new times, this seems to be a mistake or placeholder in the output.
    - `*_BigO`: Coefficient for the high-order term in the running time.
    - `*_RMS`: Normalized root-mean square error of string comparison.

