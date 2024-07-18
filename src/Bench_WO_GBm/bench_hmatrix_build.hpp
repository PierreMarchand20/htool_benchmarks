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

class FT_Generator {
  public:
    std::vector<double> p1_permuted, p2_permuted;
    std::shared_ptr<const Cluster<htool::underlying_type<double>>> m_target_root_cluster, m_source_root_cluster;
    std::unique_ptr<GeneratorTestDoubleSymmetric> generator;

    void SetUp(int nr, int nc, char Symmetry, char UPLO, htool::underlying_type<double> epsilon, double eta) {

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