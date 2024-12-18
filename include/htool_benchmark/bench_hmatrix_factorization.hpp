#ifndef HTOOL_BENCHMARK_HMATRIX_FACTORIZATION_HPP
#define HTOOL_BENCHMARK_HMATRIX_FACTORIZATION_HPP

#include "NEW_tree_builder.hpp" // où la fonction que l'on teste se situe
// #include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/linalg/interface.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
// #include <htool/matrix/linalg/interface.hpp>
#include "hmatrix_fixture.hpp"
#include <htool/testing/generate_test_case.hpp>
#include <htool/testing/generator_input.hpp>
#include <htool/testing/generator_test.hpp>

#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std;
using namespace htool;

template <typename FixtureHMatrix>
void bench_hmatrix_factorization(char symmetry_type) {
    // declare variables
    std::vector<std::string> List_algo_type;
    std::vector<double> List_epsilon;
    std::vector<int> List_pbl_size;

    // custom parameters
    const int number_of_repetitions = 9;
    List_algo_type                  = {"Classic", "TaskBased"};
    List_epsilon                    = {1e-10, 1e-8, 1e-6, 1e-4};
    List_pbl_size                   = {1 << 11, 1 << 12, 1 << 13, 1 << 14, 1 << 15};
    double eta                      = 100;
    char trans                      = 'N'; // arg of lu_solve

    // header csv file
    std::ofstream savefile;
    savefile.open("bench_hmatrix_factorization_vs_pbl_size.csv");
    savefile << "epsilon, dim_pbl, number_of_threads, algo_type, id_rep, compression_ratio, space_saving, factorization_time (s), solve_time (s) \n";

    // cout parameters
    std::cout << " ++++++++++++++++++ Test case: ++++++++++++++++++ " << std::endl;
    std::cout << "Test_case: " << "pbl_size" << std::endl;
    std::cout << "List_algo_type: " << List_algo_type << std::endl;
    std::cout << "List_epsilon: " << List_epsilon << std::endl;
    std::cout << "List_pbl_size: " << List_pbl_size << std::endl;
    std::cout << "List_thread: " << "1" << std::endl;
    std::cout << "Number_of_repetitions: " << number_of_repetitions << std::endl;
    std::cout << std::endl;

    // computation
    for (double epsilon : List_epsilon) {
        for (int dim_pbl : List_pbl_size) {
            double List_factorization_duration[number_of_repetitions] = {0};
            double List_solving_duration[number_of_repetitions]       = {0};
            double list_compression_ratio[number_of_repetitions]      = {0};
            double list_space_saving[number_of_repetitions]           = {0};
            Matrix<double> Y_dense(dim_pbl, 1, 1);
            FixtureHMatrix fixture;
            // fixture.setup_benchmark(dim_pbl, epsilon, eta, symmetry_type);

            for (string algo_type : List_algo_type) {
                double compression_ratio = 0;
                double space_saving      = 0;

                // Setup
                fixture.setup_benchmark(dim_pbl, epsilon, eta, symmetry_type, algo_type);
                // if (algo_type == "Classic") {
                //     fixture.setup_benchmark_classic(dim_pbl, epsilon, eta, symmetry_type);
                // } else if (algo_type == "TaskBased") {
                //     fixture.setup_benchmark_task_based(dim_pbl, epsilon, eta, symmetry_type);
                // }

                // Compression ratio and space saving
                auto hmatrix_information = get_hmatrix_information(*fixture.root_hmatrix);
                compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                space_saving             = std::stod(hmatrix_information["Space_saving"]);
                // print_hmatrix_information(*fixture.root_hmatrix, std::cout);

                for (int id_rep = 0; id_rep < number_of_repetitions; id_rep++) {
                    std::chrono::steady_clock::time_point start, end;

                    if (algo_type == "Classic" || algo_type == "TaskBased") {
                        // Timer
                        start = std::chrono::steady_clock::now();
                        if (symmetry_type == 'N') {
                            lu_factorization(*fixture.root_hmatrix);
                        } else if (symmetry_type == 'S') {
                            cholesky_factorization(fixture.root_hmatrix->get_UPLO(), *fixture.root_hmatrix);
                        }
                        end = std::chrono::steady_clock::now();

                        std::chrono::duration<double> duration_facto = end - start;

                        start = std::chrono::steady_clock::now();
                        if (symmetry_type == 'N') {
                            lu_solve(trans, *fixture.root_hmatrix, Y_dense);
                        } else if (symmetry_type == 'S') {
                            cholesky_solve(fixture.root_hmatrix->get_UPLO(), *fixture.root_hmatrix, Y_dense);
                        }
                        end = std::chrono::steady_clock::now();

                        std::chrono::duration<double> duration_solve = end - start;

                        // data saving
                        savefile << epsilon << ", " << dim_pbl << ", " << "1" << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration_facto.count() << ", " << duration_solve.count() << "\n";

                        List_factorization_duration[id_rep] = duration_facto.count();
                        List_solving_duration[id_rep]       = duration_solve.count();
                        list_compression_ratio[id_rep]      = compression_ratio;
                        list_space_saving[id_rep]           = space_saving;
                    }
                    // else if (algo_type == "Dense") { // Todo : Update
                    //     // Densification
                    //     Matrix<double> HA_dense(Fixture.A->get_target_cluster().get_size(), Fixture.A->get_source_cluster().get_size());
                    //     copy_to_dense(*Fixture.A, HA_dense.data());

                    //     // Timer
                    //     start = std::chrono::steady_clock::now();
                    //     lu_factorization(HA_dense);
                    //     end                                          = std::chrono::steady_clock::now();
                    //     std::chrono::duration<double> duration_facto = end - start;

                    //     start = std::chrono::steady_clock::now();
                    //     lu_solve(trans, HA_dense, *Fixture.B_dense);
                    //     end                                          = std::chrono::steady_clock::now();
                    //     std::chrono::duration<double> duration_solve = end - start;
                    //     // data saving
                    //     savefile << epsilon << ", " << dim_pbl << ", " << "1" << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration_facto.count() << ", " << duration_solve.count() << "\n";

                    //     List_factorization_duration[id_rep] = duration_facto.count();
                    //     List_solving_duration[id_rep]       = duration_solve.count();
                    //     list_compression_ratio[id_rep]      = compression_ratio;
                    //     list_space_saving[id_rep]           = space_saving;
                    // }
                }
                // mean and stddev saving
                double mean_facto, std_dev_facto, mean_solve, std_dev_solve, mean_comp_ratio, std_dev_comp_ratio, mean_space_saving, std_dev_space_saving;

                compute_standard_deviation(List_factorization_duration, number_of_repetitions, mean_facto, std_dev_facto);
                compute_standard_deviation(List_solving_duration, number_of_repetitions, mean_solve, std_dev_solve);
                compute_standard_deviation(list_compression_ratio, number_of_repetitions, mean_comp_ratio, std_dev_comp_ratio);
                compute_standard_deviation(list_space_saving, number_of_repetitions, mean_space_saving, std_dev_space_saving);

                savefile << epsilon << ", " << dim_pbl << ", " << "1" << ", " << algo_type << ", " << "mean" << ", " << mean_comp_ratio << ", " << mean_space_saving << ", " << mean_facto << ", " << mean_solve << "\n";
                savefile << epsilon << ", " << dim_pbl << ", " << "1" << ", " << algo_type << ", " << "std_dev" << ", " << std_dev_comp_ratio << ", " << std_dev_space_saving << ", " << std_dev_facto << ", " << std_dev_solve << "\n";
            }
        }
    }
    savefile.close();
}
#endif
