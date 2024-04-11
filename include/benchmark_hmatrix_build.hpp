
#include <htool/clustering/clustering.hpp>
#include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_distributed_output.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/testing/dense_blocks_generator_test.hpp>
#include <htool/testing/generator_input.hpp>
#include <htool/testing/generator_test.hpp>
#include <htool/testing/geometry.hpp>
#include <htool/testing/partition.hpp>
#include <mpi.h>
#include <task_based_tree_builder.hpp>

using namespace std;
using namespace htool;

template <typename T, typename GeneratorTestType>
bool benchmark_hmatrix_build(int nr, int nc, bool use_local_cluster, char Symmetry, char UPLO, htool::underlying_type<T> epsilon) {

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
    ClusterTreeBuilder<htool::underlying_type<T>> recursive_build_strategy;
    // recursive_build_strategy.set_partition(partition);
    // recursive_build_strategy.set_minclustersize(2);

    std::shared_ptr<const Cluster<htool::underlying_type<T>>> source_root_cluster;
    std::shared_ptr<const Cluster<htool::underlying_type<T>>> target_root_cluster = make_shared<const Cluster<htool::underlying_type<T>>>(recursive_build_strategy.create_cluster_tree(nr, 3, p1.data(), 2, sizeWorld, partition.data()));

    if (Symmetry == 'N' && nr != nc) {
        // Geometry
        double z2 = 1 + 0.1;
        create_disk(3, z2, nc, p2.data());

        // partition
        test_partition(3, nc, p2, sizeWorld, partition);

        // Clustering
        // source_recursive_build_strategy.set_minclustersize(2);

        source_root_cluster = make_shared<const Cluster<htool::underlying_type<T>>>(recursive_build_strategy.create_cluster_tree(nc, 3, p2.data(), 2, sizeWorld, partition.data()));
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
    GeneratorTestType generator(3, nr, nc, p1_permuted, p2_permuted, *target_root_cluster, *source_root_cluster, false, false);

    // HMatrix
    double eta = 10;

    std::unique_ptr<HMatrixTaskBasedTreeBuilder<T, htool::underlying_type<T>>> hmatrix_tree_builder;
    if (use_local_cluster) {
        hmatrix_tree_builder = std::make_unique<HMatrixTaskBasedTreeBuilder<T, htool::underlying_type<T>>>(target_root_cluster->get_cluster_on_partition(rankWorld), source_root_cluster->get_cluster_on_partition(rankWorld), epsilon, eta, Symmetry, UPLO, -1, -1);
    } else {
        hmatrix_tree_builder = std::make_unique<HMatrixTaskBasedTreeBuilder<T, htool::underlying_type<T>>>(*target_root_cluster, *source_root_cluster, epsilon, eta, Symmetry, UPLO, -1, rankWorld);
    }

    // build
    auto root_hmatrix = hmatrix_tree_builder->build(generator);

    save_leaves_with_rank(root_hmatrix, "leaves_" + htool::NbrToStr(rankWorld));
    save_levels(root_hmatrix, "level_" + htool::NbrToStr(rankWorld) + "_", {0, 1, 2});

    if (rankWorld == 0) {
        print_tree_parameters(root_hmatrix, std::cout);
        print_hmatrix_information(root_hmatrix, std::cout);
    }
    print_distributed_hmatrix_information(root_hmatrix, std::cout, MPI_COMM_WORLD);

    return is_error;
}
