#include "benchmark/benchmark.h"
#include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/linalg/interface.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
#include <htool/matrix/linalg/interface.hpp>
#include <htool/testing/generate_test_case.hpp>
#include <htool/testing/generator_input.hpp>

#include <htool/testing/generator_test.hpp>

using namespace std;
using namespace htool;
using namespace benchmark::internal;
using namespace benchmark;

/* Parameters of the benchmarks */
const int min_number_of_rows = 1 << 10; // min_number_of_rows must be a power of 2 else "->Ranges(...)" will run unexpected benchmarks. Reminder : (x << y) <=> x * pow(2, y)
const int max_number_of_rows = 1 << 10; // max_number_of_rows must be a power of 2. Benchmarks will loop from (min_number_of_rows) to (max_number_of_rows) with (2) as common ratio.

const int min_number_of_threads = 1;
const int max_number_of_threads = 1;

const int min_exp_epsilon = 3;
const int max_exp_epsilon = 3;

/* Fixture */
class FT_Facto : public ::benchmark::Fixture {
  public:
    std::unique_ptr<TestCaseSolve<double, GeneratorTestDoubleSymmetric>> test_case;
    std::unique_ptr<HMatrix<double, htool::underlying_type<double>>> A;

    void SetUp(const ::benchmark::State &state) override {
        int n1         = state.range(0);
        int n2         = n1;
        double epsilon = std::pow(10, -state.range(1));
        double eta     = 100;

        // Setup test case
        test_case = std::make_unique<htool::TestCaseSolve<double, GeneratorTestDoubleSymmetric>>('L', 'N', n1, n2, 1, -1);

        // HMatrix
        HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder_A(*test_case->root_cluster_A_output, *test_case->root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
        A = std::make_unique<HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder_A.build(*test_case->operator_A));

        // Matrix
        int ni_A = test_case->root_cluster_A_input->get_size();
        int no_A = test_case->root_cluster_A_output->get_size();
        int ni_X = test_case->root_cluster_X_input->get_size();
        int no_X = test_case->root_cluster_X_output->get_size();
        Matrix<double> B_dense(no_X, ni_X);
        generate_random_matrix(B_dense);
    }
};

/* functions to benchmark */
BENCHMARK_DEFINE_F(FT_Facto, BM_Classic) // Classic implementation
(benchmark::State &state) {

    for (auto _ : state) { /*Timed zone*/
        // LU factorization
        lu_factorization(*A);
    }

    auto count = static_cast<size_t>(state.range(0)); // square matrix version
    state.SetComplexityN(count);
}

BENCHMARK_REGISTER_F(FT_Facto, BM_Classic)
    ->RangeMultiplier(2)
    ->Ranges({{min_number_of_rows, max_number_of_rows}, {min_exp_epsilon, max_exp_epsilon}}) // square matrix version
    ->ArgNames({"N", "ExpEps"})
    ->ThreadRange(min_number_of_threads, max_number_of_threads)
    ->Complexity(benchmark::oNLogN);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    MPI_Finalize();
    return 0;
}