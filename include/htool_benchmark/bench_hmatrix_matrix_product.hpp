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
    List_epsilon                    = {1e-10, 1e-8, 1e-6, 1e-4};

    if (test_case_type == "pbl_size") { // 1<<19 vs 1 thread OK sur Cholesky, 1<<20 vs 1 thread out of memory
        List_pbl_size = {1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19};
        List_thread   = {1};
    }
    if (test_case_type == "thread") {
        List_pbl_size = {1 << 19};
        List_thread   = {1, 2, 4, 8, 16};
    }
    if (test_case_type == "ratio") {
        List_pbl_size = {1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19};
        List_thread   = {1, 2, 4, 8, 16};
    }

    // header csv file
    std::ofstream savefile;
    savefile.open("bench_hmatrix_matrix_product_vs_" + test_case_type + ".csv");
    savefile << "epsilon, dim_pbl, number_of_threads, algo_type, id_rep, compression_ratio, space_saving, time (s) | mean time (s) | standard_deviation \n";

    // cout parameters
    std::cout << " ++++++++++++++++++ Test case: ++++++++++++++++++ " << std::endl;
    std::cout << "Number_of_repetitions: " << number_of_repetitions << std::endl;
    std::cout << "Test_case: " << test_case_type << std::endl;
    std::cout << "List_algo_type: " << List_algo_type << std::endl;
    std::cout << "List_epsilon: " << List_epsilon << std::endl;
    std::cout << "List_pbl_size: " << List_pbl_size << std::endl;
    std::cout << "List_thread: " << List_thread << std::endl;
    std::cout << std::endl;

    // computation
    for (double epsilon : List_epsilon) {
        id_pbl_size = 0;

        for (int dim_pbl : List_pbl_size) {
            // Setup
            FixtureHMatrix fixture;
            double eta = 10;
            if (symmetry_type != 'N') {
                fixture.setup_benchmark(dim_pbl, epsilon, eta, symmetry_type);
            } else {
                fixture.setup_benchmark(dim_pbl, dim_pbl, epsilon, eta);
            }

            double list_matrix_product_duration[number_of_repetitions] = {0};

            for (std::string algo_type : List_algo_type) {
                // max_dim_pbl <= (1 << 15) for dense on laptop else "Abandon (core dumped)"
                id_thread     = 0;
                is_ratio_done = false;

                for (int n_threads : List_thread) {
                    // To avoid crossed terms in ratio_pbl_size_thread case
                    if (test_case_type == "ratio_pbl_size_thread") {
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
                        double compression_ratio = 0;
                        double space_saving      = 0;
                        std::vector<double> x(dim_pbl), y(dim_pbl, 1);

                        if (algo_type == "Classic") {
                            // Timer
                            start = std::chrono::steady_clock::now();
                            for (int i = 0; i < number_of_products; i++) {
                                openmp_internal_add_hmatrix_vector_product('N', 1., *fixture.root_hmatrix, y.data(), 0., x.data());
                            }
                            end                                    = std::chrono::steady_clock::now();
                            std::chrono::duration<double> duration = end - start;

                            // Compression ratio and space saving
                            auto hmatrix_information = get_hmatrix_information(*fixture.root_hmatrix);
                            compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                            space_saving             = std::stod(hmatrix_information["Space_saving"]);

                            // data saving
                            savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration.count() << "\n";
                            list_matrix_product_duration[id_rep] = duration.count();

                        } else if (algo_type == "TaskBased") {
                            // Timer
                            start = std::chrono::steady_clock::now();
                            for (int i = 0; i < number_of_products; i++) {
                                NEW_openmp_add_hmatrix_vector_product('N', 1., *fixture.root_hmatrix, y.data(), 0., x.data());
                            }
                            end                                    = std::chrono::steady_clock::now();
                            std::chrono::duration<double> duration = end - start;

                            // Compression ratio and space saving
                            auto hmatrix_information = get_hmatrix_information(*fixture.root_hmatrix);
                            compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                            space_saving             = std::stod(hmatrix_information["Space_saving"]);

                            // data saving
                            savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration.count() << "\n";
                            list_matrix_product_duration[id_rep] = duration.count();
                        }
                    }
                    // mean and stddev saving
                    double mean, std_dev;
                    compute_standard_deviation(list_matrix_product_duration, number_of_repetitions, mean, std_dev);
                    savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << "mean" << ", " << "N.A" << ", " << "N.A" << ", " << mean << "\n";
                    savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << "stddev" << ", " << "N.A" << ", " << "N.A" << ", " << std_dev << "\n";

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
