#ifndef HTOOL_BENCHMARK_HMATRIX_FIXTURE_HPP
#define HTOOL_BENCHMARK_HMATRIX_FIXTURE_HPP

#include "generator_fixture.hpp"
// #include "task_based_tree_builder.hpp" // où la fonction que l'on teste se situe
#include <htool/hmatrix/tree_builder/tree_builder.hpp>

using namespace std;
using namespace htool;
template <typename FixtureGenerator>
class FixtureHMatrix {

    std::unique_ptr<FixtureGenerator> fixture_generator;

  public:
    std::unique_ptr<htool::HMatrix<double, htool::underlying_type<double>>> root_hmatrix;
    int max_node_L0 = 64;

    void setup_benchmark_classic(int N, htool::underlying_type<double> epsilon, double eta, char symmetry) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(N);

        std::unique_ptr<HMatrixTreeBuilder<double, htool::underlying_type<double>>> hmatrix_tree_builder;
        hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilder<double, htool::underlying_type<double>>>(epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L');
        root_hmatrix         = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder->build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_target_root_cluster), -1, -1, false, max_node_L0));
    }

    void setup_benchmark_classic(int M, int N, htool::underlying_type<double> epsilon, double eta) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(M, N);

        std::unique_ptr<HMatrixTreeBuilder<double, htool::underlying_type<double>>> hmatrix_tree_builder;
        hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilder<double, htool::underlying_type<double>>>(epsilon, eta, 'N', 'N');

        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder->build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_source_root_cluster), -1, -1, false, max_node_L0));
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

        std::unique_ptr<HMatrixTreeBuilder<double, htool::underlying_type<double>>> hmatrix_tree_builder;
        hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilder<double, htool::underlying_type<double>>>(epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L');

        if (benchmark_type == "TaskBased") {
            root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder->build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_target_root_cluster), -1, -1, true, max_node_L0));

        } else if (benchmark_type == "Classic") {

            root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder->build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_target_root_cluster), -1, -1, false, max_node_L0));

        } else {
            throw std::runtime_error("Unknown benchmark type");
        }
    }
};
#endif