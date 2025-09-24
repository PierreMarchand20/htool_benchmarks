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
    std::vector<HMatrix<double> *> L0;
    int max_node_L0 = 64;

    void setup_benchmark_classic(int N, htool::underlying_type<double> epsilon, double eta, char symmetry) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(N);

        HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L');
        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.openmp_build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_target_root_cluster), -1, -1));
    }

    void setup_benchmark_classic(int M, int N, htool::underlying_type<double> epsilon, double eta) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(M, N);

        HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(epsilon, eta, 'N', 'N');
        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_source_root_cluster), -1, -1));
    }

    void setup_benchmark(int N, htool::underlying_type<double> epsilon, double eta, char symmetry, std::string benchmark_type) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(N);

        HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L');

        if (benchmark_type == "TaskBased") {
            root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.openmp_build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_target_root_cluster), -1, -1));

        } else if (benchmark_type == "Classic") {
            root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.task_based_build(*(fixture_generator->generator), *(fixture_generator->m_target_root_cluster), *(fixture_generator->m_target_root_cluster), L0, max_node_L0, -1, -1));

        } else {
            throw std::runtime_error("Unknown benchmark type");
        }
    }
};
#endif
