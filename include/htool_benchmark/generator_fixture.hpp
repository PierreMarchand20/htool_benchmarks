#ifndef HTOOL_BENCHMARK_GENERATOR_FIXTURE_HPP
#define HTOOL_BENCHMARK_GENERATOR_FIXTURE_HPP

#include "htool/clustering/tree_builder/tree_builder.hpp"
// #include <memory>
// #include <vector>
// #include <htool/clustering/clustering.hpp>

namespace htool_benchmark {
template <typename GeneratorType>
class FixtureGenerator {
    std::shared_ptr<htool::ClusterTreeBuilder<double>> m_cluster_tree_builder;

  public:
    std::vector<double> p1, p2;
    std::shared_ptr<const htool::Cluster<htool::underlying_type<double>>> m_target_root_cluster, m_source_root_cluster;
    std::unique_ptr<GeneratorType> generator;

    using CoefficientPrecision = typename GeneratorType::CoefficientPrecision;

    FixtureGenerator() : m_cluster_tree_builder(std::make_shared<htool::ClusterTreeBuilder<double>>()) {}
    FixtureGenerator(std::shared_ptr<htool::ClusterTreeBuilder<double>> cluster_tree_builder) : m_cluster_tree_builder(cluster_tree_builder) {}

    void setup_benchmark(int N) {
        srand(1);

        // Geometry
        p1.resize(3 * N);
        double radius          = 1.;
        double step            = 1.75 * M_PI * radius / sqrt((double)N);
        double length          = 2 * M_PI * radius;
        double pointsPerCircle = length / step;
        double angleStep       = 2 * M_PI / pointsPerCircle;
        for (int j = 0; j < N; j++) {
            p1[3 * j + 0] = radius * cos(angleStep * j);
            p1[3 * j + 1] = radius * sin(angleStep * j);
            p1[3 * j + 2] = (step * j) / pointsPerCircle;
        }

        m_target_root_cluster = std::make_shared<const htool::Cluster<htool::underlying_type<double>>>(m_cluster_tree_builder->create_cluster_tree(N, 3, p1.data(), 2, 2));

        // Generator
        if constexpr (GeneratorType::require_permuted_input) {
            auto &permutation = m_target_root_cluster->get_permutation();
            std::vector<double> permuted_p1(p1.size());
            for (int j = 0; j < N; j++) {
                permuted_p1[j + 0]     = p1[3 * permutation[j] + 0];
                permuted_p1[j + N]     = p1[3 * permutation[j] + 1];
                permuted_p1[j + 2 * N] = p1[3 * permutation[j] + 2];
            }
            p1 = permuted_p1;
        }
        generator = std::make_unique<GeneratorType>(3, p1, p1);
    }

    void setup_benchmark(int nr, int nc) {
        srand(1);

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

        m_target_root_cluster = std::make_shared<const htool::Cluster<htool::underlying_type<double>>>(m_cluster_tree_builder->create_cluster_tree(nr, 3, p1.data(), 2, 2));

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

        m_source_root_cluster = std::make_shared<const htool::Cluster<htool::underlying_type<double>>>(m_cluster_tree_builder->create_cluster_tree(nc, 3, p2.data(), 2, 2));

        // Generator
        if constexpr (GeneratorType::require_permuted_input) {
            auto &target_permutation = m_target_root_cluster->get_permutation();
            std::vector<double> permuted_p1(p1.size());
            for (int j = 0; j < nr; j++) {
                permuted_p1[j + 0]      = p1[3 * target_permutation[j] + 0];
                permuted_p1[j + nr]     = p1[3 * target_permutation[j] + 1];
                permuted_p1[j + 2 * nr] = p1[3 * target_permutation[j] + 2];
            }
            p1 = permuted_p1;

            auto &source_permutation = m_source_root_cluster->get_permutation();
            std::vector<double> permuted_p2(p2.size());
            for (int j = 0; j < nc; j++) {
                permuted_p2[j + 0]      = p1[3 * source_permutation[j] + 0];
                permuted_p2[j + nc]     = p1[3 * source_permutation[j] + 1];
                permuted_p2[j + 2 * nc] = p1[3 * source_permutation[j] + 2];
            }
            p2 = permuted_p2;
        }

        generator = std::make_unique<GeneratorType>(3, p1, p2);
    }
};
} // namespace htool_benchmark
#endif
