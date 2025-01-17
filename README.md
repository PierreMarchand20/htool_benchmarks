# Htool benchmarks

This repository contains some performance tests used by the library [Htool](https://pmarchand.pages.math.cnrs.fr/htool_documentation/index.html). In particular, it is used to perform benchmarks with other libraries [here](https://github.com/PierreMarchand20/HMatrix_Benchmarks).

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

This section presents the results of our benchmarks on the task-based implementation of our algorithm. We conducted these benchmarks to evaluate the performance of our algorithm and identify areas for improvement. 

The results of these benchmarks are presented in graph format.
The benchmarks were conducted on different thread configurations and problem sizes, and we measured the execution times and performance of our algorithm. 

In addition, we also provide a detailed explanation of how to use the benchmarks and customize them to suit your specific needs.

### Overview
Our aim is to test the implementation of task-based parallelism compared with a conventional parallel implementation. 
To do this, we compare the execution times given by different test cases as a function of several parameters. 

The test cases are: the assembly of Hmatrices `bench_hmatrix_build`, the product of Hmatrices and matrices `bench_hmatrix_matrix_product`, and the factorization of Hmatrices `bench_hmatrix_factorization`. 

For each test case, we vary the size of the problem `dim`, the number of threads `number_of_threads`, the tolerance which controls the relative error on block approximation `epsilon` and the type of implementation `algo_type`. In addition, we repeat each benchmark several times `id_rep` to ensure that there is little spurious noise in the results. 

These are the Hmatrix compression ratio `compression_ratio` and the storage saving ratio `space_saving`, as well as one or more execution times depending on the test case. 
We save them in the form of csv files bearing the name of the test case. These files can then be read by a python script `plot_bench.py`, which takes care of the visualization.

### Usage
````bash
# Go into build directory
cd htool_benchmarks
mkdir build
cd build

# Compile one, several or all benchmarks
cmake ..
make bench_hmatrix_build
make bench_hmatrix_matrix_product
make bench_hmatrix_factorization
make build-benchmarks

# Run a test case
./bench_hmatrix_build <test_case>={0|1|2} <symmetry_type>[default=N]={N|S}
./bench_hmatrix_matrix_product <test_case>={0|1|2} <symmetry_type>[default=N]={N|S}
./bench_hmatrix_factorization <symmetry_type>[default=N]={N|S}

# View one or all csv files. These are automatically copied to their own directory in /ouput. Several execution options are available, depending on the test case. 
python ../src/plot_bench.py <csv_name> <is_log_scale>[default=False]={True|False}
../scripts/run_plot_bench.sh
````

The following options are available:
- <test_case>: specifies the test case to run. The available test cases are :
    - 0 : runs the benchmark with regard to the size of the problem,
    - 1 : runs the benchmark with regard to the number of thread,
    - 2 : runs the benchmark with regard to the ratio problem size to number of thread.
- <symmetry_type>: specifies the symmetry type of the matrix. The available symmetry types are:
    - N : no symmetry (by default value),
    - S : symmetric matrix.
- <is_log_scale>: specifies whether to use a logarithmic scale for the x-axis. The available options are :
    - True : uses a logarithmic scale,
    - False : uses a linear scale (by default value).

For example, to run the benchmark with test case 1, symmetry type S, execute the following command:
````bash
./bench_hmatrix_build 1 S
````

Note that the benchmark will output the results to a CSV file named <bench_type>_vs_<test_case>.csv in the /build directory. Then, the Python script will create a corresponding directory in /output and copy the CSV file into it. Finally, the graph result will also be saved in this directory.



### Benchmark customization
The benchmark parameters are accessible in the hpp files `bench_hmatrix_build`, `bench_hmatrix_matrix_product` and `bench_hmatrix_factorization` located in the `include/htool_benchmark` folder. The user will find the following parameters in the `custom parameters` section of the code : 
- `number_of_repetitions`: the number of repetitions of the benchmark,
- `List_algo_type` : the set of implementations,  
- `List_epsilon` : the set of epsilon values (the tolerance which controls the relative error on block approximation),  
- `eta`: eligibility constant in the admissibility condition (section 'Hierarchical clustering' [here](https://pmarchand.pages.math.cnrs.fr/htool_documentation/introduction/overview.html)),  
- `List_pbl_size`: the set of problem size values,  
- `List_thread`: the set of values for the number of threads,
- `number_of_products` : the number of Hmatrix products timed,
- `number_of_solves` : the number of linear system solves timed,
- `trans` : form of the system of equations A * X = B or A**T * X = B.

The plots of the results can also be customized by modifying the `custom_parameters` function in the `plot_bench.py` file located in the `src` folder. The user will find the following parameters : 
- `SubList_dim`: the subset of problem size values to display,  
- `SubList_thread`: the subset of values for the number of threads to display,
- `SubList_algo_type` : the subset of implementations to display,  
- `SubList_epsilon` : the subset of epsilon values to display,  
- `log_exponent` : the exponent 'a' of the logarithm for graph rescaling by N log<sup>a</sup> N in the case of Hmatrix building.
- `data` : the variable to be plotted. Note that in the factorization benchmark, the user can modify `data` to plot either the factorization time or the solve time.

### Results with htool_benchmark version `TODO`
The following results figures can be found in `output_reference/version_TODO` directory along corresponding CSV files.

#### Hmatrix building time
Custom parameters :
- `number_of_repetitions` = 9
- `List_algo_type` = {"Classic", "TaskBased"}
- `List_epsilon` = {1e-10, 1e-8, 1e-4};
- `eta` = 10;
- `List_pbl_size` = {2<sup>15</sup>, 2<sup>16</sup>, 2<sup>17</sup>, 2<sup>18</sup>, 2<sup>19</sup>};
- `List_thread` = {1, 2, 4, 8, 16};
- <symmetry_type> = N

![Hmatrix_build_vs_pbl_size](output_reference/Hmatrix_build_vs_pbl_size/32768__65536_131072_262144_524288/mean_time_vs_pbl_size.png "Hmatrix building time vs problem size with 1 thread")


![Hmatrix_build_vs_thread](output_reference/Hmatrix_build_vs_thread/1__2__4__8_16/mean_time_vs_thread.png "Hmatrix building time vs number of thread with problem size equals  2¹⁹")

![Hmatrix_build_vs_ratio](output_reference/Hmatrix_build_vs_ratio/_32768__65536_131072_262144_524288__1__2__4__8_16/mean_time_vs_ratio.png "Hmatrix building time vs ratio problem size on number of thread")

#### Hmatrix matrix product time

Custom parameters :
- `number_of_repetitions` = 9
- `number_of_products` = 30
- `List_algo_type` = {"Classic", "TaskBased"}
- `List_epsilon` = {1e-10, 1e-8, 1e-4};
- `eta` = 10;
- `List_pbl_size` = {2<sup>15</sup>, 2<sup>16</sup>, 2<sup>17</sup>, 2<sup>18</sup>, 2<sup>19</sup>};
- `List_thread` = {1, 2, 4, 8, 16};
- <symmetry_type> = N
  
![Hmatrix_matrix_product_vs_pbl_size](output_reference/Hmatrix_matrix_product_vs_pbl_size/32768__65536_131072_262144_524288/mean_time_vs_pbl_size.png "Hmatrix matrix product time vs problem size with 1 thread")

![Hmatrix_matrix_product_vs_thread](output_reference/Hmatrix_matrix_product_vs_thread/1__2__4__8_16/mean_time_vs_thread.png "Hmatrix matrix product time vs number of thread with problem size equals  2¹⁹")

![Hmatrix_matrix_product_vs_ratio](output_reference/Hmatrix_matrix_product_vs_ratio/_32768__65536_131072_262144_524288__1__2__4__8_16/mean_time_vs_ratio.png "Hmatrix matrix product time vs ratio problem size on number of thread")

#### Hmatrix factorizations LU and Cholesky

Custom parameters :
- `number_of_repetitions` = 2
- `number_of_solves` = 30
- `List_algo_type` = {"Classic", "TaskBased"}
- `List_epsilon` = {1e-10, 1e-7, 1e-4};
- `eta` = 100;
- `List_pbl_size` = {2<sup>15</sup>, 2<sup>16</sup>, 2<sup>17</sup>, 2<sup>18</sup>, 2<sup>19</sup>};
- `List_thread` = {1};
- <symmetry_type> = N
- `trans` = N

![Hmatrix_factorization_LU_vs_pbl_size](output_reference/Hmatrix_factorization_vs_pbl_size/16384__32768__65536_131072_262144/mean_time_facto_LU_vs_pbl_size.png "Hmatrix factorization LU time vs problem size with 1 thread")

![Hmatrix_factorization_Cho_vs_pbl_size](output_reference/Hmatrix_factorization_vs_pbl_size/16384__32768__65536_131072_262144/mean_time_facto_Cho_vs_pbl_size.png "Hmatrix factorization Cholesky time vs problem size with 1 thread")



# TODO
- Développer le contexte : expliquer pourquoi vous avez effectué ces benchmarks et quels sont les objectifs que vous souhaitez atteindre.
- Inclure des graphiques :
  - expliquer comment les lire et les reproduire.
  - analyser les résultats.
  - parler du stddev.
- mettre un numéro de version pour htool_benchmark et update ce readme avec.