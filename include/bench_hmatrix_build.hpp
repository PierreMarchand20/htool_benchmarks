#ifndef HTOOL_BENCHMARKS_HMATRIX_BUILD_HPP
#define HTOOL_BENCHMARKS_HMATRIX_BUILD_HPP

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

class FT_Generator {
  public:
    std::vector<double> p1, p2;
    std::shared_ptr<const htool::Cluster<htool::underlying_type<double>>> m_target_root_cluster, m_source_root_cluster;
    std::unique_ptr<htool::GeneratorTestDoubleSymmetric> generator;

    void SetUp(int nr, int nc, char Symmetry, char UPLO, htool::underlying_type<double> epsilon, double eta) {
        srand(1);

        // Geometry
        double z1 = 1;
        p1.resize(3 * nr);
        p2.resize(Symmetry == 'N' ? 3 * nc : 1);
        htool::create_disk(3, z1, nr, p1.data());

        // Clustering
        htool::ClusterTreeBuilder<htool::underlying_type<double>> recursive_build_strategy;

        m_target_root_cluster = std::make_shared<const htool::Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nr, 3, p1.data(), 2, 2));

        if (Symmetry == 'N' && nr != nc) {
            // Geometry
            double z2 = 1 + 0.1;
            htool::create_disk(3, z2, nc, p2.data());

            // Clustering
            // source_recursive_build_strategy.set_minclustersize(2);

            m_source_root_cluster = std::make_shared<const htool::Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nc, 3, p2.data(), 2, 2));
        } else {
            m_source_root_cluster = m_target_root_cluster;
            p2                    = p1;
        }

        // Generator
        generator = std::make_unique<htool::GeneratorTestDoubleSymmetric>(3, p1, p2);
    }

    void SetUpBench0(int nr, int nc, char Symmetry, char UPLO, htool::underlying_type<double> epsilon, double eta) {
        srand(1);

        // Clustering
        htool::ClusterTreeBuilder<htool::underlying_type<double>> recursive_build_strategy;

        // Geometry
        p1.resize(3 * nr);
        double radius          = 1.;
        double step            = 1.75 * M_PI * radius / sqrt((double)nr);
        double length          = 2 * M_PI * radius;
        double pointsPerCircle = length / step;
        double angleStep       = 2 * M_PI / pointsPerCircle;
        for (int j = 0; j < nr; j++) {
            p1[3 * j + 0] = radius * cos(angleStep * j);
            p1[3 * j + 1] = radius * sin(angleStep * j);
            p1[3 * j + 2] = (step * j) / pointsPerCircle;
        }

        m_target_root_cluster = std::make_shared<const htool::Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nr, 3, p1.data(), 2, 2));

        p2.resize(3 * nc);
        double shift    = 10;
        step            = 1.75 * M_PI * radius / sqrt((double)nc);
        pointsPerCircle = length / step;
        angleStep       = 2 * M_PI / pointsPerCircle;
        for (int j = 0; j < nc; j++) {
            p2[3 * j + 0] = shift + radius * cos(angleStep * j);
            p2[3 * j + 1] = shift + radius * sin(angleStep * j);
            p2[3 * j + 2] = shift + (step * j) / pointsPerCircle;
        }

        m_source_root_cluster = std::make_shared<const htool::Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nc, 3, p2.data(), 2, 2));

        // Generator
        generator = std::make_unique<htool::GeneratorTestDoubleSymmetric>(3, p1, p2);
    }
};

#endif