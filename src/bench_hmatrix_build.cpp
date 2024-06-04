/*
TODO :
- résoudre pbl de variance qd on boucle sur la taille du pbl
- faire plot
- Tester différente complexité, best guest so far : NLogN ?
- ajouter un setup (et un teardown)
- Màj ReadMe
- deux types de BM :
    - temps/mémoire en fonction du nombre de threads, en prenant l'exemple de la comparaison pour la
      boucle d'assemblage parallélisé avec les taches ou avec guided,
    - temps/mémoire en fonction de la taille du problème en comparant matrice dense vs hmatrix.
- ajouter la boucle suivante: for (auto epsilon : {1e-14, 1e-6}) {...}
- continuer à optimiser le build_hmatrix avec le task-based
- profiler build_hmatrix avec au moins deux profilers de differents types
-----------------------------------------------------------------------------------------------------------
*/

/*
Tuto:
- Profiling : https://medium.com/distributed-knowledge/optimizations-for-c-multi-threaded-programs-33284dee5e9c
- Doc Gbenchmark : https://github.com/google/benchmark/blob/main/docs/user_guide.md#output-formats
- Random interleaving to reduce run-to-run variance : https://github.com/google/benchmark/blob/main/docs/random_interleaving.md
  and problematic : https://github.com/google/benchmark/issues/1051
- Benchmarking tips : https://llvm.org/docs/Benchmarking.html
-------------------------------------------------------------------------------------------------------------
*/

#include "../external/htool/tests/functional_tests/hmatrix/test_hmatrix_build.hpp"
#include "NEW_hmatrix_build.hpp"
#include "benchmark/benchmark.h"

using namespace std;
using namespace htool;

/* In_Arguments of the benchmarks, must be a power of 2 else "->Ranges" will run unexpected benchmarks with a power of 2 in In_Arguments */
const int number_of_rows    = 128;
const int number_of_columns = 128;
// const int number_of_rows_increased    = 256;
// const int number_of_columns_increased = 256;
const int Pow              = 0;   // benchmarks will loop from (In_Arguments) to (2^Pow*In_Arguments) with (RangeMutiplier) as common ratio. Reminder : (x << y) <=> x * pow(2, y)
const double MinWarmUpTime = 0.2; // warmup time in order to decrease variance
const int Threads          = 1;   // number of threads
const int Repetitions      = 9;   // number of repetions per benchmark != Iterations. Must be >= 9 to perform U-test.
const double MinTime       = 0.1; // per-repetition time

/* function to benchmark */
static void BM_test_hmatrix_build(benchmark::State &state) { // Classic implementation

    // auto count = static_cast<size_t>(state.range(0) * state.range(1)); // rectangular matrix version
    auto count = static_cast<size_t>(state.range(0)); // square matrix version

    for (auto _ : state) {

        /*Setup generator and hmatrix tree builder. Not timed zone*/
        state.PauseTiming(); // Stop timers.

        // Parameters
        int nr                                 = state.range(0);
        int nc                                 = state.range(0); // square matrix version
        bool use_local_cluster                 = true;
        char Symmetry                          = 'N';
        char UPLO                              = 'N';
        htool::underlying_type<double> epsilon = 1e-14;
        bool use_dense_blocks_generator        = true;
        VirtualGenerator<double> *generator;
        // HMatrixTreeBuilder<double, htool::underlying_type<double>> *hmatrix_tree_builder;

        // Get the number of processes
        int sizeWorld;
        MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

        // Get the rankWorld of the process
        int rankWorld;
        MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);

        srand(1);
        bool is_error = false;

        // Geometry
        double z1 = 1;
        vector<double> p1(3 * nr), p1_permuted, off_diagonal_p1;
        vector<double> p2(Symmetry == 'N' ? 3 * nc : 1), p2_permuted, off_diagonal_p2;
        create_disk(3, z1, nr, p1.data());

        // Partition
        std::vector<int> partition{};
        test_partition(3, nr, p1, sizeWorld, partition);

        // Clustering
        ClusterTreeBuilder<htool::underlying_type<double>> recursive_build_strategy;
        // recursive_build_strategy.set_partition(partition);
        // recursive_build_strategy.set_minclustersize(2);

        std::shared_ptr<const Cluster<htool::underlying_type<double>>> source_root_cluster;
        std::shared_ptr<const Cluster<htool::underlying_type<double>>> target_root_cluster = make_shared<const Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nr, 3, p1.data(), 2, sizeWorld, partition.data()));

        if (Symmetry == 'N' && nr != nc) {
            // Geometry
            double z2 = 1 + 0.1;
            create_disk(3, z2, nc, p2.data());

            // partition
            test_partition(3, nc, p2, sizeWorld, partition);

            // Clustering
            // source_recursive_build_strategy.set_minclustersize(2);

            source_root_cluster = make_shared<const Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nc, 3, p2.data(), 2, sizeWorld, partition.data()));
        } else {
            source_root_cluster = target_root_cluster;
            p2                  = p1;
        }

        // Permutation on geometry
        p1_permuted.resize(3 * nr);
        const auto &target_permutation = target_root_cluster->get_permutation();
        for (int i = 0; i < target_permutation.size(); i++) {
            p1_permuted[i * 3 + 0] = p1[target_permutation[i] * 3 + 0];
            p1_permuted[i * 3 + 1] = p1[target_permutation[i] * 3 + 1];
            p1_permuted[i * 3 + 2] = p1[target_permutation[i] * 3 + 2];
        }
        p2_permuted.resize(3 * nc);
        if (Symmetry == 'N' && nr != nc) {
            const auto &source_permutation = source_root_cluster->get_permutation();
            for (int i = 0; i < source_permutation.size(); i++) {
                p2_permuted[i * 3 + 0] = p2[source_permutation[i] * 3 + 0];
                p2_permuted[i * 3 + 1] = p2[source_permutation[i] * 3 + 1];
                p2_permuted[i * 3 + 2] = p2[source_permutation[i] * 3 + 2];
            }
        } else {
            p2_permuted = p1_permuted;
        }

        // Generator
        generator = new GeneratorTestDoubleSymmetric(3, nr, nc, p1_permuted, p2_permuted, *target_root_cluster, *source_root_cluster, false, false);

        // HMatrix
        double eta = 10;

        std::unique_ptr<HMatrixTreeBuilder<double, htool::underlying_type<double>>> hmatrix_tree_builder;
        if (use_local_cluster) {
            hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilder<double, htool::underlying_type<double>>>(target_root_cluster->get_cluster_on_partition(rankWorld), source_root_cluster->get_cluster_on_partition(rankWorld), epsilon, eta, Symmetry, UPLO, -1, -1, -1);
        } else {
            hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilder<double, htool::underlying_type<double>>>(*target_root_cluster, *source_root_cluster, epsilon, eta, Symmetry, UPLO, -1, rankWorld, rankWorld);
        }

        std::shared_ptr<VirtualDenseBlocksGenerator<double>> dense_blocks_generator;
        if (use_dense_blocks_generator) {
            dense_blocks_generator = std::make_shared<DenseBlocksGeneratorTest<double>>(*generator);
        }
        hmatrix_tree_builder->set_dense_blocks_generator(dense_blocks_generator);

        state.ResumeTiming(); // And resume timers.

        /*Timed zone*/
        auto root_hmatrix = hmatrix_tree_builder->build(*generator);
        // test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(1), true, 'N', 'N', 1e-14, true); // rectangular matrix version
        // test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(0), true, 'N', 'N', 1e-14, true); // square matrix version
        // test_hmatrix_build_bis(state.range(7), state.range(8));

        /*Delete generator and hmatrix_tree_builder. Not timed zone*/
        state.PauseTiming(); // Stop timers.
        delete generator;
        state.ResumeTiming();
    };
    state.SetComplexityN(count);
    state.SetItemsProcessed(count * state.iterations());
    state.SetBytesProcessed(count * state.iterations() * sizeof(int32_t));
}
BENCHMARK(BM_test_hmatrix_build)
    ->RangeMultiplier(2)
    // ->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}, {number_of_columns, (1 << Pow) * number_of_columns}}) // rectangular matrix version
    ->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}}) // square matrix version
    // ->Arg(number_of_rows)       // 0
    // ->Arg(number_of_columns)    // 1
    // ->Arg(true)                 // 2
    // ->Arg('N')                  // 3
    // ->Arg('N')                  // 4
    // ->Arg(1e-14)                // 5
    // ->Arg(true)                 // 6
    // ->Arg(generator)            // 7
    // ->Arg(hmatrix_tree_builder) // 8
    ->MinWarmUpTime(MinWarmUpTime)
    ->Repetitions(Repetitions)
    ->Threads(Threads)
    ->Complexity(benchmark::oNLogN)
    ->MinTime(MinTime)
    // ->Setup(DoSetupOld)
    ;

/* Define benchmark of the new implementation */
static void BM_NEW_test_hmatrix_build(benchmark::State &state) { // task based implementation
    // auto count = static_cast<size_t>(state.range(0) * state.range(1)); // rectangular matrix version
    auto count = static_cast<size_t>(state.range(0)); // square matrix version

    for (auto _ : state) {
        // NEW_test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(1), true, 'N', 'N', 1e-14, true);
        NEW_test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(state.range(0), state.range(0), true, 'N', 'N', 1e-14, true); // square matrix version
    };
    state.SetComplexityN(count);
    state.SetItemsProcessed(count * state.iterations());
    state.SetBytesProcessed(count * state.iterations() * sizeof(int32_t));
}
BENCHMARK(BM_NEW_test_hmatrix_build)
    ->RangeMultiplier(2)
    // ->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}, {number_of_columns, (1 << Pow) * number_of_columns}}) // rectangular matrix version
    ->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}}) // square matrix version
    ->MinWarmUpTime(MinWarmUpTime)
    ->Repetitions(Repetitions)
    ->Threads(Threads)
    ->Complexity(benchmark::oNLogN)
    ->MinTime(MinTime);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    MPI_Finalize();
    return 0;
}