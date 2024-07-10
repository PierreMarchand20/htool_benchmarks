
#include "bench_hmatrix_matrix_product.hpp"

/* Parameters of the benchmarks */
const int min_number_of_rows = 1 << 10; // min_number_of_rows must be a power of 2 else "->Ranges(...)" will run unexpected benchmarks. Reminder : (x << y) <=> x * pow(2, y)
const int max_number_of_rows = 1 << 13; // max_number_of_rows must be a power of 2. Benchmarks will loop from (min_number_of_rows) to (max_number_of_rows) with (2) as common ratio.

const int min_number_of_threads = 1;
const int max_number_of_threads = min_number_of_threads;

const int min_exp_epsilon = 4;
const int max_exp_epsilon = 4;

BENCHMARK_REGISTER_F(FT_LinearAlgebra, BM_Classic)
    ->RangeMultiplier(2)
    ->Ranges({{min_number_of_rows, max_number_of_rows}, {min_exp_epsilon, max_exp_epsilon}}) // square matrix version
    ->ArgNames({"N", "ExpEps"})
    ->ThreadRange(min_number_of_threads, max_number_of_threads)
    ->Complexity(benchmark::oNLogN);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    MPI_Finalize();
    return 0;
}