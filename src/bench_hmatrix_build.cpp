// TODO :
// - Màj ReadMe
// - réduire au mieux la variance https://github.com/google/benchmark/blob/main/docs/reducing_variance.md
// - ajouter la boucle suivante: for (auto epsilon : {1e-14, 1e-6}) {...}
// - deux types de BM :
//     - temps/mémoire en fonction du nombre de threads, en prenant l'exemple de la comparaison pour la
//       boucle d'assemblage parallélisé avec les taches ou avec guided,
//     - temps/mémoire en fonction de la taille du problème en comparant matrice dense vs hmatrix.
// - lancer 10 fois chaque BM et faire la moyenne des runtimes car le calcul en parallel fait varier les temps
// - profiler build_hmatrix avec au moins deux profilers de differents types
// - continuer à optimiser le build_hmatrix avec le task-based

// Tuto:
// Profiling : https://medium.com/distributed-knowledge/optimizations-for-c-multi-threaded-programs-33284dee5e9c
// Doc Gbenchmark : https://github.com/google/benchmark/blob/main/docs/user_guide.md#output-formats

#include "../external/htool/tests/functional_tests/hmatrix/test_hmatrix_build.hpp"
#include "NEW_hmatrix_build.hpp"
#include "benchmark/benchmark.h"

using namespace std;
using namespace htool;

// In_Arguments of the benchmarks, must be a power of 2 else "->Ranges" will run unexpected benchmarks with a power of 2 in In_Arguments
const int number_of_rows              = 128;
const int number_of_rows_increased    = 256;
const int number_of_columns           = 128;
const int number_of_columns_increased = 256;

const int Pow = 1; // benchmarks will loop from (In_Arguments) to (2^Pow*In_Arguments) with (2) as common ratio

// Define Benchmarks

// Simple benchmark for BM_test_hmatrix_build with unique call : {
// static void BM_NEW_test_hmatrix_build(benchmark::State &state) {
//     for (auto _ : state) {
//         NEW_test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, true, 'N', 'N', 1e-14);
//     }
// }
// BENCHMARK(BM_NEW_test_hmatrix_build);

// static void BM_test_hmatrix_build(benchmark::State &state) {
//     for (auto _ : state) {
//         test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, true, 'N', 'N', 1e-14, true);
//     }
// }
// BENCHMARK(BM_test_hmatrix_build);
// }

// More complex benchmark for BM_NEW_test_hmatrix_build and BM_test_hmatrix_build with multiple calls : {
// static void BM_NEW_test_hmatrix_build(benchmark::State &state) {
//     for (auto _ : state) {
//         NEW_test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(1), true, 'N', 'N', 1e-14);
//     }
// }
// // BENCHMARK(BM_NEW_test_hmatrix_build)->Args({number_of_rows, number_of_columns}); // Equivalent to previous "Simple benchmark for BM_test_hmatrix_build with unique call"
// BENCHMARK(BM_NEW_test_hmatrix_build)->RangeMultiplier(2)->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}, {number_of_columns, (1 << Pow) * number_of_columns}}); // loop from (In_Arguments) to (2^Pow*In_Arguments) with (2) as common ratio

static void BM_test_hmatrix_build(benchmark::State &state) {
    for (auto _ : state) {
        test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(1), true, 'N', 'N', 1e-14, true);
    }
}
BENCHMARK(BM_test_hmatrix_build)->RangeMultiplier(2)->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}, {number_of_columns, (1 << Pow) * number_of_columns}});
// }

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();

    MPI_Finalize();
    return 0;
}