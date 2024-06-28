/*
TODO :
- résoudre pbl de variance qd on boucle sur la taille du pbl
- faire plot
- Tester différente complexité, best guest so far : NLogN ?
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
- Google Benchmark basic guide : https://ccfd.github.io/courses/hpc_lab01.html#Learning_Google_benchmark
-------------------------------------------------------------------------------------------------------------
*/

// #include "../external/htool/tests/functional_tests/hmatrix/test_hmatrix_build.hpp"
// #include "NEW_hmatrix_build.hpp"
#include "benchmark/benchmark.h"

#include "NEW_tree_builder.hpp" // où la fonction que l'on teste se situe
#include <htool/clustering/clustering.hpp>
#include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_distributed_output.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp> // où la fonction que l'on teste se situe
#include <htool/testing/dense_blocks_generator_test.hpp>
#include <htool/testing/generator_input.hpp>
#include <htool/testing/generator_test.hpp>
#include <htool/testing/geometry.hpp>
#include <htool/testing/partition.hpp>
#include <mpi.h>

using namespace std;
using namespace htool;
using namespace benchmark::internal;
using namespace benchmark;

/* In_Arguments of the benchmarks, must be a power of 2 else "->Ranges" will run unexpected benchmarks with a power of 2 in In_Arguments */
const int number_of_rows = 128;
// const int number_of_columns = 128;
// const int number_of_rows_increased    = 256;
// const int number_of_columns_increased = 256;

const int Pow     = 1; // benchmarks will loop from (In_Arguments) to (2^Pow*In_Arguments) with (RangeMutiplier) as common ratio. Reminder : (x << y) <=> x * pow(2, y)
const int Threads = 1; // number of threads

/* Fixture */
class FT_Generator : public ::benchmark::Fixture {
  public:
    std::vector<double> p1_permuted, p2_permuted;
    std::shared_ptr<const Cluster<htool::underlying_type<double>>> m_target_root_cluster, m_source_root_cluster;
    std::unique_ptr<GeneratorTestDoubleSymmetric> generator;

    void SetUp(const ::benchmark::State &state) override {
        // Parameters
        int nr                                 = state.range(0);
        int nc                                 = state.range(0); // square matrix version
        char Symmetry                          = 'N';
        char UPLO                              = 'N';
        htool::underlying_type<double> epsilon = 1e-14;
        double eta                             = 10;

        // Get the number of processes
        int sizeWorld;
        MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

        // Get the rankWorld of the process
        int rankWorld;
        MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);

        srand(1);
        // bool is_error = false;

        // Geometry
        double z1 = 1;
        vector<double> p1(3 * nr);
        vector<double> p2(Symmetry == 'N' ? 3 * nc : 1);
        create_disk(3, z1, nr, p1.data());

        // Partition
        std::vector<int> partition{};
        test_partition(3, nr, p1, sizeWorld, partition);

        // Clustering
        ClusterTreeBuilder<htool::underlying_type<double>> recursive_build_strategy;
        // recursive_build_strategy.set_partition(partition);
        // recursive_build_strategy.set_minclustersize(2);

        m_target_root_cluster = make_shared<const Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nr, 3, p1.data(), 2, sizeWorld, partition.data()));

        if (Symmetry == 'N' && nr != nc) {
            // Geometry
            double z2 = 1 + 0.1;
            create_disk(3, z2, nc, p2.data());

            // partition
            test_partition(3, nc, p2, sizeWorld, partition);

            // Clustering
            // source_recursive_build_strategy.set_minclustersize(2);

            m_source_root_cluster = make_shared<const Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nc, 3, p2.data(), 2, sizeWorld, partition.data()));
        } else {
            m_source_root_cluster = m_target_root_cluster;
            p2                    = p1;
        }

        // Permutation on geometry
        p1_permuted.resize(3 * nr);
        const auto &target_permutation = m_target_root_cluster->get_permutation();
        for (int i = 0; i < target_permutation.size(); i++) {
            p1_permuted[i * 3 + 0] = p1[target_permutation[i] * 3 + 0];
            p1_permuted[i * 3 + 1] = p1[target_permutation[i] * 3 + 1];
            p1_permuted[i * 3 + 2] = p1[target_permutation[i] * 3 + 2];
        }
        p2_permuted.resize(3 * nc);
        if (Symmetry == 'N' && nr != nc) {
            const auto &source_permutation = m_source_root_cluster->get_permutation();
            for (int i = 0; i < source_permutation.size(); i++) {
                p2_permuted[i * 3 + 0] = p2[source_permutation[i] * 3 + 0];
                p2_permuted[i * 3 + 1] = p2[source_permutation[i] * 3 + 1];
                p2_permuted[i * 3 + 2] = p2[source_permutation[i] * 3 + 2];
            }
        } else {
            p2_permuted = p1_permuted;
        }

        // Generator
        generator = std::make_unique<GeneratorTestDoubleSymmetric>(3, nr, nc, p1_permuted, p2_permuted, *m_target_root_cluster, *m_source_root_cluster, false, false);
    }
};

/* functions to benchmark */
BENCHMARK_DEFINE_F(FT_Generator, BM_Classic) // Classic implementation
(benchmark::State &state) {

    char Symmetry                          = 'N';
    char UPLO                              = 'N';
    htool::underlying_type<double> epsilon = 1e-14;
    double eta                             = 10;

    using HMatrixTreeBuilderType = HMatrixTreeBuilder<double, htool::underlying_type<double>>;

    for (auto _ : state) { /*Timed zone*/
        // Hmatrix
        std::unique_ptr<HMatrixTreeBuilderType> hmatrix_tree_builder;
        hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilderType>(m_target_root_cluster->get_cluster_on_partition(0), m_source_root_cluster->get_cluster_on_partition(0), epsilon, eta, Symmetry, UPLO, -1, -1, -1);

        auto root_hmatrix = hmatrix_tree_builder->build(*generator);
    }

    auto count = static_cast<size_t>(state.range(0)); // square matrix version
    state.SetComplexityN(count);
}

BENCHMARK_REGISTER_F(FT_Generator, BM_Classic)
    ->RangeMultiplier(2)
    ->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}}) // square matrix version
    ->ArgName({"N"})
    ->Threads(Threads)
    ->Complexity(benchmark::oNLogN);

BENCHMARK_DEFINE_F(FT_Generator, BM_TaskBased) // Task based implementation
(benchmark::State &state) {

    char Symmetry                          = 'N';
    char UPLO                              = 'N';
    htool::underlying_type<double> epsilon = 1e-14;
    double eta                             = 10;

    using HMatrixTreeBuilderType = HMatrixTaskBasedTreeBuilder<double, htool::underlying_type<double>>;

    for (auto _ : state) { /*Timed zone*/
        // Hmatrix
        std::unique_ptr<HMatrixTreeBuilderType> hmatrix_tree_builder;
        hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilderType>(m_target_root_cluster->get_cluster_on_partition(0), m_source_root_cluster->get_cluster_on_partition(0), epsilon, eta, Symmetry, UPLO, -1, -1, -1);

        auto root_hmatrix = hmatrix_tree_builder->build(*generator);
    }

    auto count = static_cast<size_t>(state.range(0)); // square matrix version
    state.SetComplexityN(count);
}

BENCHMARK_REGISTER_F(FT_Generator, BM_TaskBased)
    ->RangeMultiplier(2)
    ->Ranges({{number_of_rows, (1 << Pow) * number_of_rows}}) // square matrix version
    ->ArgName({"N"})
    ->Threads(Threads)
    ->Complexity(benchmark::oNLogN);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    MPI_Finalize();
    return 0;
}