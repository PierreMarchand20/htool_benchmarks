#ifndef HTOOL_BENCHMARK_HMATRIX_FIXTURE_HPP
#define HTOOL_BENCHMARK_HMATRIX_FIXTURE_HPP

#include "generator_fixture.hpp"
#include "task_based_tree_builder.hpp" // où la fonction que l'on teste se situe
#include <htool/hmatrix/tree_builder/tree_builder.hpp>

template <typename FixtureGenerator>
class FixtureHMatrix {

    std::unique_ptr<FixtureGenerator> fixture_generator;

  public:
    std::unique_ptr<htool::HMatrix<double, htool::underlying_type<double>>> root_hmatrix;

    void setup_benchmark_classic(int N, htool::underlying_type<double> epsilon, double eta, char symmetry) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(N);
        htool::HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*fixture_generator->m_target_root_cluster, *fixture_generator->m_target_root_cluster, epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L', -1, -1, -1);
        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*(fixture_generator->generator)));
    }

    void setup_benchmark_classic(int M, int N, htool::underlying_type<double> epsilon, double eta) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(M, N);
        htool::HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*fixture_generator->m_target_root_cluster, *fixture_generator->m_source_root_cluster, epsilon, eta, 'N', 'N', -1, -1, -1);
        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*(fixture_generator->generator)));
    }

    // void setup_benchmark_task_based(int N, htool::underlying_type<double> epsilon, double eta, char symmetry) {
    //     fixture_generator = std::make_unique<FixtureGenerator>();
    //     fixture_generator->setup_benchmark(N);
    //     htool::HMatrixTaskBasedTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*fixture_generator->m_target_root_cluster, *fixture_generator->m_target_root_cluster, epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L', -1, -1, -1);
    //     root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*(fixture_generator->generator)));
    // }

    void setup_benchmark(int N, htool::underlying_type<double> epsilon, double eta, char symmetry, std::string benchmark_type) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(N);
        if (benchmark_type == "TaskBased") {
            htool::HMatrixTaskBasedTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*fixture_generator->m_target_root_cluster, *fixture_generator->m_target_root_cluster, epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L', -1, -1, -1);
            root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*(fixture_generator->generator)));

        } else if (benchmark_type == "Classic") {
            htool::HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*fixture_generator->m_target_root_cluster, *fixture_generator->m_target_root_cluster, epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L', -1, -1, -1);
            root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*(fixture_generator->generator)));
        } else {
            throw std::runtime_error("Unknown benchmark type");
        }
    }
};
#endif