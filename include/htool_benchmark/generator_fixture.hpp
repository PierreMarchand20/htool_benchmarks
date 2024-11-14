#ifndef HTOOL_BENCHMARK_GENERATOR_FIXTURE_HPP
#define HTOOL_BENCHMARK_GENERATOR_FIXTURE_HPP

#include <htool/clustering/cluster.hpp>

using namespace std;   // TODO: a enlever
using namespace htool; // TODO: a enlever

namespace htool_benchmark {
template <typename GeneratorType>
class FixtureGenerator {
  public:
    std::vector<double> p1, p2;
    std::shared_ptr<const Cluster<htool::underlying_type<double>>> m_target_root_cluster, m_source_root_cluster;
    std::unique_ptr<GeneratorType> generator;

    void setup_benchmark(int N) {
        srand(1);

        // Clustering
        ClusterTreeBuilder<htool::underlying_type<double>> recursive_build_strategy;

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

        m_target_root_cluster = make_shared<const Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(N, 3, p1.data(), 2, 2));

        // Generator
        generator = std::make_unique<GeneratorType>(3, p1, p1);
    }

    void setup_benchmark(int nr, int nc) {
        srand(1);

        // Clustering
        ClusterTreeBuilder<htool::underlying_type<double>> recursive_build_strategy;

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

        m_target_root_cluster = make_shared<const Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nr, 3, p1.data(), 2, 2));

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

        m_source_root_cluster = make_shared<const Cluster<htool::underlying_type<double>>>(recursive_build_strategy.create_cluster_tree(nc, 3, p2.data(), 2, 2));

        // Generator
        generator = std::make_unique<GeneratorType>(3, p1, p2);
    }
};
} // namespace htool_benchmark
#endif