#ifndef HTOOL_BENCHMARKS_HMATRIX_MATRIX_PRODUCT_HPP
#define HTOOL_BENCHMARKS_HMATRIX_MATRIX_PRODUCT_HPP

#include "NEW_add_hmatrix_vector_product.hpp" // où la fonction que l'on teste se situe
#include <htool/hmatrix/linalg/add_hmatrix_vector_product.hpp>

#include "hmatrix_fixture.hpp"
#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
#include <iostream>

using namespace htool;

namespace htool_benchmark {

/**
 * @brief Benchmark H-matrix product.
 *
 * This function benchmarks the product of an H-matrix with a vector using
 * different algorithms, problem sizes, and precision levels.
 *
 * @param test_case_type Type of test case: "pbl_size", "thread" or "ratio".
 * @param symmetry_type Type of symmetry for the H-matrix: 'N' for non-symmetric,
 * 'S' for symmetric, and 'H' for hermitian.
 */
template <typename FixtureHMatrix>
void bench_hmatrix_matrix_product(std::string test_case_type, char symmetry_type) {

    // declare variables
    std::vector<std::string> List_algo_type;
    std::vector<double> List_epsilon;
    std::vector<int> List_pbl_size;
    std::vector<int> List_thread;
    int id_pbl_size(0), id_thread(0);
    bool is_ratio_done(false);

    // custom parameters
    const int number_of_repetitions = 9;
    const int number_of_products    = 30;
    List_algo_type                  = {"Classic", "TaskBased"};
    List_epsilon                    = {1e-10, 1e-7, 1e-4};
    double eta                      = 10;

    if (test_case_type == "pbl_size") { // 1<<19 vs 1 thread OK sur Cholesky, 1<<20 vs 1 thread out of memory
        List_pbl_size = {1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19};
        // List_pbl_size = {1 << 10};
        List_thread = {1};
    }
    if (test_case_type == "thread") {
        List_pbl_size = {1 << 19};
        List_thread   = {1, 2, 4, 8, 16};
    }
    if (test_case_type == "ratio") {
        List_pbl_size = {1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19};
        // List_pbl_size = {1 << 10};
        List_thread = {1, 2, 4, 8, 16};
    }

    // header csv file
    std::ofstream savefile;
    savefile.open("bench_hmatrix_matrix_product_vs_" + test_case_type + ".csv");
    savefile << "epsilon, dim, number_of_threads, algo_type, id_rep, compression_ratio, space_saving, time (s) \n";

    // cout parameters
    std::cout << " ++++++++++++++++++ Test case: ++++++++++++++++++ " << std::endl;
    std::cout << "Test_case: " << test_case_type << std::endl;
    std::cout << "List_algo_type: " << List_algo_type << std::endl;
    std::cout << "List_epsilon: " << List_epsilon << std::endl;
    std::cout << "List_pbl_size: " << List_pbl_size << std::endl;
    std::cout << "List_thread: " << List_thread << std::endl;
    std::cout << "Number_of_repetitions: " << number_of_repetitions << std::endl;
    std::cout << "Number_of_products: " << number_of_products << std::endl;
    std::cout << "Symmetry_type: " << symmetry_type << std::endl;
    std::cout << "Eta: " << eta << std::endl;
    std::cout << std::endl;

    // computation
    for (double epsilon : List_epsilon) {
        id_pbl_size = 0;

        for (int dim : List_pbl_size) {
            // Setup
            FixtureHMatrix fixture;
            if (symmetry_type != 'N') {
                fixture.setup_benchmark_classic(dim, epsilon, eta, symmetry_type);
            } else {
                fixture.setup_benchmark_classic(dim, dim, epsilon, eta);
            }

            double list_matrix_product_duration[number_of_repetitions] = {0};
            double list_compression_ratio[number_of_repetitions]       = {0};
            double list_space_saving[number_of_repetitions]            = {0};

            for (std::string algo_type : List_algo_type) { // max_dim <= (1 << 15) for dense on laptop else "Abandon (core dumped)"
                id_thread     = 0;
                is_ratio_done = false;

                for (int n_threads : List_thread) {
                    // To avoid crossed terms in ratio_pbl_size_thread case
                    if (test_case_type == "ratio") {
                        if (is_ratio_done) {
                            is_ratio_done = false;
                            break;
                        }
                        if (id_pbl_size != id_thread) {
                            id_thread++;
                            continue;
                        }
                    }
                    omp_set_num_threads(n_threads);

                    for (int id_rep = 0; id_rep < number_of_repetitions; id_rep++) {
                        std::chrono::steady_clock::time_point start, end;
                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(*fixture.root_hmatrix);
                        double compression_ratio = std::stod(hmatrix_information["Compression_ratio"]);
                        double space_saving      = std::stod(hmatrix_information["Space_saving"]);
                        // print_hmatrix_information(*fixture.root_hmatrix, std::cout);
                        std::vector<double> x(dim), y(dim, 1);
                        std::chrono::duration<double> duration;
                        // Timer
                        if (algo_type == "Classic") {
                            start = std::chrono::steady_clock::now();
                            for (int i = 0; i < number_of_products; i++) {
                                openmp_internal_add_hmatrix_vector_product('N', 1., *fixture.root_hmatrix, y.data(), 0., x.data());
                            }
                            end = std::chrono::steady_clock::now();

                        } else if (algo_type == "TaskBased") {
                            start = std::chrono::steady_clock::now();
                            for (int i = 0; i < number_of_products; i++) {
                                NEW_openmp_add_hmatrix_vector_product('N', 1., *fixture.root_hmatrix, y.data(), 0., x.data());
                            }
                            end = std::chrono::steady_clock::now();
                        } else {
                            std::cerr << "Unknown algo_type: " << algo_type << std::endl;
                            exit(1);
                        }

                        duration = end - start;

                        // data saving
                        savefile << epsilon << ", " << dim << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration.count() << "\n";
                        list_matrix_product_duration[id_rep] = duration.count();
                        list_compression_ratio[id_rep]       = compression_ratio;
                        list_space_saving[id_rep]            = space_saving;
                    }
                    // mean and stddev saving
                    double mean_prod, std_dev_prod, mean_comp_ratio, std_dev_comp_ratio, mean_space_saving, std_dev_space_saving;

                    compute_standard_deviation(list_matrix_product_duration, number_of_repetitions, mean_prod, std_dev_prod);
                    compute_standard_deviation(list_compression_ratio, number_of_repetitions, mean_comp_ratio, std_dev_comp_ratio);
                    compute_standard_deviation(list_space_saving, number_of_repetitions, mean_space_saving, std_dev_space_saving);

                    savefile << epsilon << ", " << dim << ", " << n_threads << ", " << algo_type << ", " << "mean" << ", " << mean_comp_ratio << ", " << mean_space_saving << ", " << mean_prod << "\n";
                    savefile << epsilon << ", " << dim << ", " << n_threads << ", " << algo_type << ", " << "stddev" << ", " << std_dev_comp_ratio << ", " << std_dev_space_saving << ", " << std_dev_prod << "\n";

                    is_ratio_done = true;
                }
            }
            id_pbl_size++;
        }
    }
    savefile.close();
}
}; // namespace htool_benchmark
#endif
