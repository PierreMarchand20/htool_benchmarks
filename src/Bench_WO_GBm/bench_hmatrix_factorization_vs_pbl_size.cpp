
#include "NEW_tree_builder.hpp" // où la fonction que l'on teste se situe
#include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/linalg/interface.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
#include <htool/matrix/linalg/interface.hpp>
#include <htool/testing/generate_test_case.hpp>
#include <htool/testing/generator_input.hpp>
#include <htool/testing/generator_test.hpp>

#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std;
using namespace htool;

/* Fixture */
class FT_Facto {
  public:
    std::unique_ptr<TestCaseSolve<double, GeneratorTestDoubleSymmetric>> test_case;
    std::unique_ptr<HMatrix<double, htool::underlying_type<double>>> A;

    void SetUp(int n1, int n2, htool::underlying_type<double> epsilon, double eta) {
        // Setup test case
        test_case = std::make_unique<htool::TestCaseSolve<double, GeneratorTestDoubleSymmetric>>('L', 'N', n1, n2, 1, -1);

        // // HMatrix
        // HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder_A(*test_case->root_cluster_A_output, *test_case->root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
        // A = std::make_unique<HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder_A.build(*test_case->operator_A));

        // // Matrix
        // int ni_A = test_case->root_cluster_A_input->get_size();
        // int no_A = test_case->root_cluster_A_output->get_size();
        // int ni_X = test_case->root_cluster_X_input->get_size();
        // int no_X = test_case->root_cluster_X_output->get_size();
        // Matrix<double> B_dense(no_X, ni_X);
        // generate_random_matrix(B_dense);
    }
};

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    const int number_of_repetitions = 9;
    const int min_dim_pbl           = 1 << 10;
    const int max_dim_pbl           = 1 << 13;
    const int dim_pbl_step          = 2;

    std::ofstream savefile;
    savefile.open("bench_hmatrix_factorization_vs_pbl_size.csv");
    savefile << "epsilon, dim_pbl, algo_type, id_rep, time (s) | mean time (s) | standard_deviation \n";

    for (double epsilon : {1e-14, 1e-10, 1e-6}) {
        for (int dim_pbl = min_dim_pbl; dim_pbl <= max_dim_pbl; dim_pbl *= dim_pbl_step) {
            // Setup
            FT_Facto Fixture;
            double eta = 100;
            Fixture.SetUp(dim_pbl, dim_pbl, epsilon, eta);
            double List_factorization_duration[number_of_repetitions] = {0};

            for (string algo_type : {"Dense", "Classic", "TaskBased"}) {
                for (int id_rep = 0; id_rep < number_of_repetitions; id_rep++) {
                    std::chrono::steady_clock::time_point start, end;
                    if (algo_type == "Classic" || algo_type == "Dense") { // HMatrixTreeBuilder
                        // HMatrix
                        HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder_A(*Fixture.test_case->root_cluster_A_output, *Fixture.test_case->root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
                        Fixture.A = std::make_unique<HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder_A.build(*Fixture.test_case->operator_A));

                        // Timer
                        start = std::chrono::steady_clock::now();
                        lu_factorization(*Fixture.A);
                        end                                    = std::chrono::steady_clock::now();
                        std::chrono::duration<double> duration = end - start;

                        // data saving
                        savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << id_rep << ", " << duration.count() << "\n";
                        List_factorization_duration[id_rep] = duration.count();

                    } else { // HMatrixTaskBasedTreeBuilder
                        // HMatrix
                        HMatrixTaskBasedTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder_A(*Fixture.test_case->root_cluster_A_output, *Fixture.test_case->root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
                        Fixture.A = std::make_unique<HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder_A.build(*Fixture.test_case->operator_A));

                        // Timer
                        start = std::chrono::steady_clock::now();
                        lu_factorization(*Fixture.A);
                        end                                    = std::chrono::steady_clock::now();
                        std::chrono::duration<double> duration = end - start;

                        // data saving
                        savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << id_rep << ", " << duration.count() << "\n";
                        List_factorization_duration[id_rep] = duration.count();
                    }
                }
                // mean and stddev saving
                double mean, std_dev;
                compute_standard_deviation(List_factorization_duration, number_of_repetitions, mean, std_dev);
                savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << "mean" << ", " << mean << "\n";
                savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << "stddev" << ", " << std_dev << "\n";
            }
        }
    }
    savefile.close();
    MPI_Finalize();
    return 0;
}