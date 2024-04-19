// TODO :
// 1) pouvoir lancer les deux benchmarks avec/sans task based en même temps. Actuellemnt redifinition de fonction qd il y a les deux include .hpp
// 2) ajouter au script bash le changement de mode du CPU pour corriger : "***WARNING*** CPU scaling is enabled, the benchmark real time measurements may be noisy and will incur extra overhead."
// 3) ajouter les boucles suivantes:
//      for (auto use_local_cluster : {true, false}) {
//          for (auto epsilon : {1e-14, 1e-6}) {
//              for (auto use_dense_block_generator : {true, false}) {

// #include "../external/htool/tests/functional_tests/hmatrix/test_hmatrix_build.hpp"
#include "NEW_hmatrix_build.hpp"
#include "benchmark/benchmark.h"

using namespace std;
using namespace htool;

const int number_of_rows              = 200;
const int number_of_rows_increased    = 400;
const int number_of_columns           = 200;
const int number_of_columns_increased = 400;

void BM_benchmark_hmatrix_build(benchmark::State &state) {
    for (auto _ : state) {
        NEW_test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, true, 'N', 'N', 1e-14);
    }
}
BENCHMARK(BM_benchmark_hmatrix_build);

// void BM_test_hmatrix_build(benchmark::State &state) {
//     for (auto _ : state) {
//         test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, true, 'N', 'N', 1e-14, true);
//     }
// }
// BENCHMARK(BM_test_hmatrix_build);

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    ::benchmark::Initialize(&argc, argv);

    ::benchmark::RunSpecifiedBenchmarks();

    MPI_Finalize();
    return 0;
}