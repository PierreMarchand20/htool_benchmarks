
#include "benchmark/benchmark.h"

#include <htool/clustering/clustering.hpp>
#include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/linalg/interface.hpp>
#include <htool/hmatrix/lrmat/SVD.hpp>
#include <htool/hmatrix/lrmat/linalg/interface.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
#include <htool/testing/generate_test_case.hpp>
#include <htool/testing/generator_input.hpp>
#include <htool/testing/generator_test.hpp>
#include <htool/testing/geometry.hpp>
#include <htool/testing/partition.hpp>

using namespace std;
using namespace htool;
using namespace benchmark::internal;
using namespace benchmark;

/* Fixture */
class FT_LinearAlgebra : public ::benchmark::Fixture {
  public:
    std::unique_ptr<TestCaseProduct<double, GeneratorTestDoubleSymmetric>> test_case;
    std::vector<double> B_vec, C_vec;
    double alpha, beta;
    char transa = 'N';
    std::unique_ptr<HMatrix<double, htool::underlying_type<double>>> root_hmatrix;

    void SetUp(const ::benchmark::State &state) override {
        double epsilon = std::pow(10, -state.range(1));
        double eta     = 10;
        test_case      = std::make_unique<TestCaseProduct<double, GeneratorTestDoubleSymmetric>>(transa, 'N', state.range(0), state.range(0), 1, 1, 2, -1);

        const Cluster<htool::underlying_type<double>> *root_cluster_A_output, *root_cluster_A_input, *root_cluster_B_output, *root_cluster_B_input, *root_cluster_C_output, *root_cluster_C_input;

        root_cluster_A_output = &test_case->root_cluster_A_output->get_cluster_on_partition(0);
        root_cluster_A_input  = &test_case->root_cluster_A_input->get_cluster_on_partition(0);
        root_cluster_B_output = &test_case->root_cluster_B_output->get_cluster_on_partition(0);
        root_cluster_B_input  = &test_case->root_cluster_B_input->get_cluster_on_partition(0);
        root_cluster_C_output = &test_case->root_cluster_C_output->get_cluster_on_partition(0);
        root_cluster_C_input  = &test_case->root_cluster_C_input->get_cluster_on_partition(0);

        HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*root_cluster_A_output, *root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
        hmatrix_tree_builder.set_low_rank_generator(std::make_shared<SVD<double>>());

        // build
        root_hmatrix = std::make_unique<HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*test_case->operator_A));

        // Dense matrix
        int ni_A = root_cluster_A_input->get_size();
        int no_A = root_cluster_A_output->get_size();
        int ni_B = root_cluster_B_input->get_size();
        int no_B = root_cluster_B_output->get_size();
        int ni_C = root_cluster_C_input->get_size();
        int no_C = root_cluster_C_output->get_size();

        // Random Input matrix
        B_vec.resize(no_B);
        C_vec.resize(no_C);
        generate_random_vector(B_vec);
        generate_random_vector(C_vec);

        generate_random_scalar(alpha);
        generate_random_scalar(beta);
    }
};

/* functions to benchmark */
BENCHMARK_DEFINE_F(FT_LinearAlgebra, BM_Classic) // Classic implementation
(benchmark::State &state) {

    for (auto _ : state) { /*Timed zone*/
        openmp_add_hmatrix_vector_product(transa, alpha, *root_hmatrix, B_vec.data(), beta, C_vec.data());
    }

    auto count = static_cast<size_t>(state.range(0)); // square matrix version
    state.SetComplexityN(count);
}
