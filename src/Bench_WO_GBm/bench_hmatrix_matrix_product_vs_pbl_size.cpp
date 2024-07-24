#include "NEW_tree_builder.hpp" // où la fonction que l'on teste se situe
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

#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std;
using namespace htool;

/* Fixture */
class FT_LinearAlgebra {
  public:
    std::unique_ptr<TestCaseProduct<double, GeneratorTestDoubleSymmetric>> test_case;
    std::vector<double> B_vec, C_vec;
    double alpha, beta;
    std::unique_ptr<HMatrix<double, htool::underlying_type<double>>> root_hmatrix, root_hmatrix_task_based;
    std::unique_ptr<Matrix<double>> HA_dense;

    void SetUp(int n1, int n2, htool::underlying_type<double> epsilon, double eta, char transa) {
        test_case = std::make_unique<TestCaseProduct<double, GeneratorTestDoubleSymmetric>>(transa, 'N', n1, n2, 1, 1, 2, -1);

        const Cluster<htool::underlying_type<double>> *root_cluster_A_output, *root_cluster_A_input, *root_cluster_B_output, *root_cluster_B_input, *root_cluster_C_output, *root_cluster_C_input;

        root_cluster_A_output = &test_case->root_cluster_A_output->get_cluster_on_partition(0);
        root_cluster_A_input  = &test_case->root_cluster_A_input->get_cluster_on_partition(0);
        root_cluster_B_output = &test_case->root_cluster_B_output->get_cluster_on_partition(0);
        // root_cluster_B_input  = &test_case->root_cluster_B_input->get_cluster_on_partition(0);
        root_cluster_C_output = &test_case->root_cluster_C_output->get_cluster_on_partition(0);
        // root_cluster_C_input  = &test_case->root_cluster_C_input->get_cluster_on_partition(0);

        // HMatrixTreeBuilder
        HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*root_cluster_A_output, *root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
        root_hmatrix = std::make_unique<HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*test_case->operator_A));

        HA_dense = std::make_unique<Matrix<double>>(root_hmatrix->get_target_cluster().get_size(), root_hmatrix->get_source_cluster().get_size());
        copy_to_dense(*root_hmatrix, HA_dense->data());

        // HMatrixTaskBasedTreeBuilder
        HMatrixTaskBasedTreeBuilder<double, htool::underlying_type<double>> hmatrix_task_based_tree_builder(*root_cluster_A_output, *root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
        root_hmatrix_task_based = std::make_unique<HMatrix<double, htool::underlying_type<double>>>(hmatrix_task_based_tree_builder.build(*test_case->operator_A));

        // Dense matrix
        int no_B = root_cluster_B_output->get_size();
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

int main(int argc, char **argv) {
    const int number_of_repetitions = 9;
    const int min_dim_pbl           = 1 << 10;
    const int max_dim_pbl           = 1 << 11;
    const int dim_pbl_step          = 2;

    std::ofstream savefile;
    savefile.open("bench_hmatrix_matrix_product_vs_pbl_size.csv");
    savefile << "epsilon, dim_pbl, algo_type, id_rep, compression_ratio, space_saving, time (s) | mean time (s) | standard_deviation \n";

    for (double epsilon : {1e-14, 1e-10}) {
        for (int dim_pbl = min_dim_pbl; dim_pbl <= max_dim_pbl; dim_pbl *= dim_pbl_step) {
            // Setup
            FT_LinearAlgebra Fixture;
            double eta  = 10;
            char transa = 'N';
            Fixture.SetUp(dim_pbl, dim_pbl, epsilon, eta, transa);
            double List_duration[number_of_repetitions] = {0};

            for (string algo_type : {"Dense", "Classic", "TaskBased"}) {
                for (int id_rep = 0; id_rep < number_of_repetitions; id_rep++) {
                    std::chrono::steady_clock::time_point start, end;
                    double compression_ratio = 0;
                    double space_saving      = 0;
                    if (algo_type == "Classic") {
                        // Timer
                        start = std::chrono::steady_clock::now();
                        openmp_add_hmatrix_vector_product(transa, Fixture.alpha, *Fixture.root_hmatrix, Fixture.B_vec.data(), Fixture.beta, Fixture.C_vec.data());
                        end                                    = std::chrono::steady_clock::now();
                        std::chrono::duration<double> duration = end - start;

                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(*Fixture.root_hmatrix);
                        compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                        space_saving             = std::stod(hmatrix_information["Space_saving"]);

                        // data saving
                        savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration.count() << "\n";
                        List_duration[id_rep] = duration.count();
                    } else if (algo_type == "TaskBased") {
                        // Timer
                        start = std::chrono::steady_clock::now();
                        openmp_add_hmatrix_vector_product(transa, Fixture.alpha, *Fixture.root_hmatrix_task_based, Fixture.B_vec.data(), Fixture.beta, Fixture.C_vec.data());
                        end                                    = std::chrono::steady_clock::now();
                        std::chrono::duration<double> duration = end - start;

                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(*Fixture.root_hmatrix_task_based);
                        compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                        space_saving             = std::stod(hmatrix_information["Space_saving"]);

                        // data saving
                        savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration.count() << "\n";
                        List_duration[id_rep] = duration.count();
                    } else if (algo_type == "Dense") {
                        // Timer
                        start = std::chrono::steady_clock::now();
                        add_matrix_vector_product(transa, Fixture.alpha, *Fixture.HA_dense, Fixture.B_vec.data(), Fixture.beta, Fixture.C_vec.data());
                        end                                    = std::chrono::steady_clock::now();
                        std::chrono::duration<double> duration = end - start;

                        // data saving
                        savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << id_rep << ", " << "N.A." << ", " << "N.A." << ", " << duration.count() << "\n";
                        List_duration[id_rep] = duration.count();
                    }
                }
                // mean and stddev saving
                double mean, std_dev;
                compute_standard_deviation(List_duration, number_of_repetitions, mean, std_dev);
                savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << "mean" << ", " << "N.A." << ", " << "N.A." << ", " << mean << "\n";
                savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << "stddev" << ", " << "N.A." << ", " << "N.A." << ", " << std_dev << "\n";
            }
        }
    }
    savefile.close();
    return 0;
}