# Htool benchmarks

This repository contains some performance tests used by the library [Htool](https://pmarchand.pages.math.cnrs.fr/htool_documentation/index.html). In particular, it is used to perform benchmarks with other libraries [here](https://github.com/PierreMarchand20/HMatrix_Benchmarks).

**ToDo Results can be find ...**

We conducted these benchmarks to evaluate the performance of our algorithm and identify areas for improvement, results are presented in graph format.
The benchmarks were conducted on different thread configurations and problem sizes, and we measured the execution times and performance of our algorithm.

In addition, we also provide a detailed explanation of how to use the benchmarks and customize them to suit your specific needs.

## Overview

Our aim is to test the Hmatrix implementation in the Htool library. To do this, we study the execution times given by different functions with several sets of parameters. The studied functions correspond to the following type of benchmark:

- The assembly of Hmatrices `bench_hmatrix_build`,
- The product of Hmatrices and matrices `bench_hmatrix_matrix_product`,
- The factorization of Hmatrices `bench_hmatrix_factorization`.

Each feature is tested against at least one of the following test cases:

- The problem size `_vs_pbl_size`,
- The number of thread `_vs_thread`,
- The ratio of the problem size on the number of thread (weak scaling) `_vs_ratio`.

The results are saved as CSV files. The above-mentioned names are used in the executables and CSV files names. For example: `bench_hmatrix_build_vs_pbl_size`, `bench_hmatrix_matrix_product_vs_pbl_size`, `bench_hmatrix_matrix_product_vs_thread`, etc...

In the result CSV files one can find the following columns:

- The tolerance which controls the relative error on block approximation `epsilon`,
- The size of the problem `dim`,
- The number of threads `number_of_threads`,
- The type of implementation `algo_type` (see remark below), 
- The index of the repetition `id_rep` of the benchmark,
- The compression ratio `compression_ratio`,
- The storage saving ratio `space_saving`,
- The execution time `time (s)`. If several times are measured, several columns are added with appropriated name, for example: `factorization_time (s)` and `solve_time (s)` in Hmatrix factorization benchmarks.
  
**Remark** : The ‘Classic’ implementation corresponds to the usual use of shared parallelism with OpenMP, typically the `# pragma omp for` instruction. The ‘TaskBased’ implementation, on the other hand, uses OpenMP's task parallelism features with the `# pragma omp taskloop` instruction.

These CSV files can then be read by a python script `plot_bench.py`, which takes care of the visualization.

## Usage

````bash
# Go into build directory
mkdir build
cd build

# Compile all benchmarks
cmake ..
make build-benchmarks

# Run a benchmark based on a specific test case
./bench_hmatrix_build <test_case>={0|1|2} <symmetry_type>[default=N]={N|S}
./bench_hmatrix_matrix_product <test_case>={0|1|2} <symmetry_type>[default=N]={N|S}
./bench_hmatrix_factorization <symmetry_type>[default=N]={N|S}

# Result CSV files are automatically copied to their own directory in ouput/ before visualization.
python ../src/plot_bench.py <csv_name> <is_log_scale>[default=False]={True|False} # View one specific CSV file
../scripts/run_plot_bench.sh # View all CSV files specified in run_plot_bench.sh
````

The following execution options are available:

- <test_case>: specifies the test case to run. The available test cases are :
  - 0 : runs the benchmark with regard to the size of the problem `_vs_pbl_size`,
  - 1 : runs the benchmark with regard to the number of thread `_vs_thread`,
  - 2 : runs the benchmark with regard to the problem size to number of thread ratio `_vs_ratio`.
- <symmetry_type>: specifies the symmetry type of the matrix. The available symmetry types are:
  - N : no symmetry (by default value),
  - S : symmetric matrix.
- <is_log_scale>: specifies whether to use a logarithmic scale for the x-axis. The available options are :
  - False : uses a linear scale (by default value),
  - True : uses a logarithmic scale.

For example, to run the `bench_hmatrix_build` benchmark with test case `1` and symmetry type `S`, execute the following command:

````bash
./bench_hmatrix_build 1 S
````

Note that the benchmark will output the results to a CSV file named `<bench_type>_vs_<test_case>.csv` in the `build/` directory. Then, the Python script will create a corresponding directory in the `output/` directory and copy the CSV file into it. Finally, the graph result will also be saved in this directory.

## Benchmark customization

The benchmark parameters are accessible in the hpp files `bench_hmatrix_build`, `bench_hmatrix_matrix_product` and `bench_hmatrix_factorization` located in the `include/htool_benchmark` folder. The user will find the following parameters in the `custom parameters` section of the code :

- `number_of_repetitions`: the number of repetitions of the benchmark (see remark below),
- `List_algo_type` : the set of implementations,  
- `List_epsilon` : the set of epsilon values (the tolerance which controls the relative error on block approximation),  
- `eta`: eligibility constant in the admissibility condition (section 'Hierarchical clustering' [here](https://pmarchand.pages.math.cnrs.fr/htool_documentation/introduction/overview.html)),  
- `List_pbl_size`: the set of problem size values,  
- `List_thread`: the set of values for the number of threads,
- `number_of_products` : the number of Hmatrix products timed,
- `number_of_solves` : the number of linear system solves timed,
- `trans` : form of the system of equations A * X = B or A**T * X = B.

**Remark**: We recommend that the number of repetitions `number_of_repetitions` is at least 9. This will enable the user to detect and correct any noise in the measured execution times. In fact, the standard deviation of the times measured over all the repetitions appears in the graphs.

The plots of the results can also be customized by modifying the `custom_parameters` function in the `plot_bench.py` file located in the `src/` folder. The user will find the following parameters :

- `SubList_dim`: the subset of problem size values to display,  
- `SubList_thread`: the subset of values for the number of threads to display,
- `SubList_algo_type` : the subset of implementations to display,  
- `SubList_epsilon` : the subset of epsilon values to display,  
- `log_exponent` : the exponent 'a' of the logarithm for graph rescaling by N log<sup>a</sup> N in the test cases `_vs_pbl_size`.
- `data` : the execution time to be plotted. Note that in the factorization benchmark, the user can modify `data` to plot either the factorization time or the solve time. 

## Results with htool_benchmark version `TODO`

The following results figures can be found in `output_reference/version_TODO` directory along corresponding CSV files.

## Hmatrix building time

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

## Hmatrix matrix product time

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

## Hmatrix factorizations LU and Cholesky

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

## TODO

- Inclure des graphiques :
  - expliquer comment les lire et les reproduire.
  - analyser les résultats.
  - parler du stddev.
- mettre un numéro de version pour htool_benchmark et update ce readme avec.