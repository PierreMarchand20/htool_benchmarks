#include "bench_hmatrix_build.hpp"

/* Parameters of the benchmarks */
const int min_number_of_rows = 1 << 10; // min_number_of_rows must be a power of 2 else "->Ranges(...)" will run unexpected benchmarks. Reminder : (x << y) <=> x * pow(2, y)
const int max_number_of_rows = 1 << 13; // max_number_of_rows must be a power of 2. Benchmarks will loop from (min_number_of_rows) to (max_number_of_rows) with (2) as common ratio.

const int min_number_of_threads = 1;
const int max_number_of_threads = 1;

BENCHMARK_REGISTER_F(FT_Generator, BM_Classic)
    ->RangeMultiplier(2)
    ->Ranges({{min_number_of_rows, max_number_of_rows}}) // square matrix version
    ->ArgName({"N"})
    ->ThreadRange(min_number_of_threads, max_number_of_threads)
    ->Complexity(benchmark::oNLogN);

BENCHMARK_REGISTER_F(FT_Generator, BM_TaskBased)
    ->RangeMultiplier(2)
    ->Ranges({{min_number_of_rows, max_number_of_rows}}) // square matrix version
    ->ArgName({"N"})
    ->ThreadRange(min_number_of_threads, max_number_of_threads)
    ->Complexity(benchmark::oNLogN);

BENCHMARK_REGISTER_F(FT_Generator, BM_Dense)
    ->RangeMultiplier(2)
    ->Ranges({{min_number_of_rows, max_number_of_rows}}) // square matrix version
    ->ArgName({"N"})
    ->ThreadRange(min_number_of_threads, max_number_of_threads)
    ->Complexity(benchmark::oNLogN);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    MPI_Finalize();
    return 0;
}