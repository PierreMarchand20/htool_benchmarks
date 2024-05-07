// TODO :
// - Màj ReadMe
// - ajouter la boucle suivante: for (auto epsilon : {1e-14, 1e-6}) {...}
// - ajouter un setup (et un teardown)
// - deux types de BM :
//     - temps/mémoire en fonction du nombre de threads, en prenant l'exemple de la comparaison pour la
//       boucle d'assemblage parallélisé avec les taches ou avec guided,
//     - temps/mémoire en fonction de la taille du problème en comparant matrice dense vs hmatrix.
// - lancer 10 fois chaque BM et faire la moyenne des runtimes car le calcul en parallel fait varier les temps
// - Tester utilité de DoNotOptimize and ClobberMemory ?
// - Tester différente complexité, best guest so far : NLogN ?
// - profiler build_hmatrix avec au moins deux profilers de differents types
// - continuer à optimiser le build_hmatrix avec le task-based

// Tuto:
// Profiling : https://medium.com/distributed-knowledge/optimizations-for-c-multi-threaded-programs-33284dee5e9c
// Doc Gbenchmark : https://github.com/google/benchmark/blob/main/docs/user_guide.md#output-formats
// Random interleaving to reduce run-to-run variance : https://github.com/google/benchmark/blob/main/docs/random_interleaving.md and problematic : https://github.com/google/benchmark/issues/1051
// Benchmarking tips : https://llvm.org/docs/Benchmarking.html

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

const int Pow              = 0;    // benchmarks will loop from (In_Arguments) to (2^Pow*In_Arguments) with (RangeMutiplier) as common ratio. Reminder : (x << y) <=> x * pow(2, y)
const double MinWarmUpTime = 0.25; // warmup time in order to decrease variance
const int Threads          = 8;    // number of threads
const int Repetitions      = 15;   // number of repetions per benchmark != Iterations.
const double MinTime       = 0.1;  // per-repetition time

// Define Benchmarks

static void BM_test_hmatrix_build(benchmark::State &state) {
    auto count = static_cast<size_t>(state.range(0) * state.range(1));
    for (auto _ : state) {
        test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(1), true, 'N', 'N', 1e-14, true);
    };
    state.SetComplexityN(count);
    state.SetItemsProcessed(count * state.iterations());
    state.SetBytesProcessed(count * state.iterations() * sizeof(int32_t));
}
BENCHMARK(BM_test_hmatrix_build)->RangeMultiplier(2)->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}, {number_of_columns, (1 << Pow) * number_of_columns}})->MinWarmUpTime(MinWarmUpTime)->Repetitions(Repetitions)->Threads(Threads)->Complexity(benchmark::oNLogN)->MinTime(MinTime);

static void BM_NEW_test_hmatrix_build(benchmark::State &state) {
    auto count = static_cast<size_t>(state.range(0) * state.range(1));
    for (auto _ : state) {
        NEW_test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(1), true, 'N', 'N', 1e-14, true);
    };
    state.SetComplexityN(count);
    state.SetItemsProcessed(count * state.iterations());
    state.SetBytesProcessed(count * state.iterations() * sizeof(int32_t));
}
BENCHMARK(BM_NEW_test_hmatrix_build)->RangeMultiplier(2)->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}, {number_of_columns, (1 << Pow) * number_of_columns}})->MinWarmUpTime(MinWarmUpTime)->Repetitions(Repetitions)->Threads(Threads)->Complexity(benchmark::oNLogN)->MinTime(MinTime);

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();

    MPI_Finalize();
    return 0;
}