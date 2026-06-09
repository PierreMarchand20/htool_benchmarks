#ifndef HTOOL_BENCHMARK_HMATRIX_FACTORIZATION_HPP
#define HTOOL_BENCHMARK_HMATRIX_FACTORIZATION_HPP

// #include "task_based_tree_builder.hpp" // où la fonction que l'on teste se situe
// #include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/linalg/factorization.hpp>
#include <htool/hmatrix/linalg/task_based_factorization.hpp> // for task_based_lu_factorization
#include <htool/hmatrix/tree_builder/tree_builder.hpp>

#include "cli.hpp"
#include "htool/hmatrix/execution_policies.hpp"
#include "htool_benchmark/utils.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <htool/testing/generate_test_case.hpp>
#include <htool/testing/generator_input.hpp>
#include <htool/testing/generator_test.hpp>
#include <iostream>

using namespace std;
using namespace htool;

/**
 * @brief Benchmark H-matrix factorization.
 *
 * This function benchmarks the factorization and solve process of H-matrices
 * using different algorithms, problem sizes, and precision levels.
 *
 * @tparam FixtureHMatrix The fixture type used for setting up H-matrix benchmarks.
 * @param symmetry_type Type of symmetry for the H-matrix: 'N' for non-symmetric, 'S' for symmetric.
 */
template <template <template <typename> class, typename> class FixtureHMatrix, template <typename> class FixtureGenerator, typename GeneratorType>
void bench_hmatrix_factorization(std::string test_case_type, char symmetry_type, std::string generator_type, std::string clustering_type, std::string low_rank_generator_type, std::string hardware_type, std::string version) {
    using CoefficientPrecision = typename FixtureGenerator<GeneratorType>::CoefficientPrecision;

    // declare variables
    std::vector<std::string> List_policy_type;
    std::vector<double> List_epsilon;
    std::vector<int> List_pbl_size;
    std::vector<int> List_thread;
    int id_pbl_size(0), id_thread(0);
    bool is_ratio_done(false);

    // custom parameters
    const int number_of_repetitions = 1;
    const int number_of_solves      = 30;
    List_policy_type                = {"seq", "par"};
    List_epsilon                    = {1e-10, 1e-7, 1e-4};
    double eta                      = 10;
    char trans                      = 'N'; // arg of lu_solve

    if (test_case_type == "pbl_size") {
        // List_pbl_size = {1 << 12, 1 << 13, 1 << 14}; // OK up to 1<<19 with Cholesky on my laptop
        // List_pbl_size = {1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19};
        List_pbl_size = {1 << 9};
        List_thread   = {1};
    }
    if (test_case_type == "thread") {
        // List_pbl_size = {1 << 19};
        // List_thread   = {1, 2, 4, 8, 16};
    }
    if (test_case_type == "ratio") {
        // List_pbl_size = {1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19};
        // // List_pbl_size = {1 << 10};
        // List_thread = {1, 2, 4, 8, 16};
    }

    // header csv file
    std::string filename = "bench_hmatrix_factorization_vs_" + test_case_type + ".csv";
    std::ofstream savefile;
    bool file_already_exists = std::filesystem::exists(filename);
    savefile.open(filename, std::ios::app);
    if (!file_already_exists) {
        savefile << "epsilon,size,number_of_threads,policy_type,id_rep,compression_ratio,space_saving, factorization_time,solve_time,clustering_type,low_rank_generator_type,symmetry_type,generator_type,hardware_type,version\n";
    }

    // cout parameters
    std::cout << " ++++++++++++++++++ Test case: ++++++++++++++++++ " << std::endl;
    std::cout << "Test_case: " << test_case_type << std::endl;
    std::cout << "List_policy_type: " << List_policy_type << std::endl;
    std::cout << "List_epsilon: " << List_epsilon << std::endl;
    std::cout << "List_pbl_size: " << List_pbl_size << std::endl;
    std::cout << "List_thread: " << List_thread << std::endl;
    std::cout << "Number_of_repetitions: " << number_of_repetitions << std::endl;
    std::cout << "Symmetry_type: " << symmetry_type << std::endl;
    std::cout << "Generator_type: " << generator_type << std::endl;
    std::cout << "Clustering_type: " << clustering_type << std::endl;
    std::cout << "Low_rank_generator_type: " << low_rank_generator_type << std::endl;
    std::cout << "Hardware_type: " << hardware_type << std::endl;
    std::cout << "Version: " << version << std::endl;
    std::cout << "Eta: " << eta << std::endl;
    std::cout << "Trans: " << trans << std::endl;
    std::cout << std::endl;

    // Partitioning strategy
    auto partitioning_strategy = process_clustering_type<double>(clustering_type);
    auto cluster_tree_builder  = std::make_shared<htool::ClusterTreeBuilder<double>>();
    cluster_tree_builder->set_partitioning_strategy(partitioning_strategy);

    // computation
    for (double epsilon : List_epsilon) {
        id_pbl_size = 0;

        for (int dim : List_pbl_size) {
            // Setup
            std::shared_ptr<FixtureGenerator<GeneratorType>> generator_fixture = std::make_shared<FixtureGenerator<GeneratorType>>(cluster_tree_builder);

            FixtureHMatrix hmatrix_fixture(generator_fixture);
            hmatrix_fixture.setup_benchmark(dim, epsilon, eta, symmetry_type, low_rank_generator_type);

            double List_factorization_duration[number_of_repetitions] = {0};
            double List_solving_duration[number_of_repetitions]       = {0};
            double list_compression_ratio[number_of_repetitions]      = {0};
            double list_space_saving[number_of_repetitions]           = {0};

            for (string policy_type : List_policy_type) {
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
                        auto local_hmatrix = *hmatrix_fixture.root_hmatrix;
                        std::chrono::steady_clock::time_point start, end;
                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(*hmatrix_fixture.root_hmatrix);
                        double compression_ratio = std::stod(hmatrix_information["Compression_ratio"]);
                        double space_saving      = std::stod(hmatrix_information["Space_saving"]);
                        // print_hmatrix_information(*fixture.root_hmatrix, std::cout);

                        Matrix<CoefficientPrecision> Y_dense(dim, 1, 1);

                        if (policy_type == "seq") {
                            start = std::chrono::steady_clock::now();
                            if (symmetry_type == 'N') {
                                lu_factorization(exec_compat::seq, local_hmatrix);
                            } else if (symmetry_type == 'S') {
                                cholesky_factorization(exec_compat::seq, local_hmatrix.get_UPLO(), local_hmatrix);
                            }
                            end = std::chrono::steady_clock::now();

                            std::chrono::duration<double> duration_facto = end - start;

                            // Timer for solve
                            start = std::chrono::steady_clock::now();
                            for (int i = 0; i < number_of_solves; i++) { // in place solve so Y_dense is variable
                                if (symmetry_type == 'N') {
                                    lu_solve(trans, local_hmatrix, Y_dense);
                                } else if (symmetry_type == 'S') {
                                    cholesky_solve(local_hmatrix.get_UPLO(), local_hmatrix, Y_dense);
                                }
                            }
                            end = std::chrono::steady_clock::now();

                            std::chrono::duration<double> duration_solve = end - start;

                            // data saving
                            savefile << epsilon << "," << dim << "," << "1" << "," << policy_type << "," << id_rep << "," << compression_ratio << "," << space_saving << "," << duration_facto.count() << "," << duration_solve.count() << "," << clustering_type << "," << low_rank_generator_type << "," << symmetry_type << "," << generator_type << "," << hardware_type << "," << version << "\n";

                            List_factorization_duration[id_rep] = duration_facto.count();
                            List_solving_duration[id_rep]       = duration_solve.count();
                            list_compression_ratio[id_rep]      = compression_ratio;
                            list_space_saving[id_rep]           = space_saving;
                        } else if (policy_type == "par" || policy_type == "omp_task") {
                            // Timer for factorization
                            start = std::chrono::steady_clock::now();
                            if (symmetry_type == 'N') {
                                lu_factorization(exec_compat::par, local_hmatrix);
                            } else if (symmetry_type == 'S') {
                                cholesky_factorization(exec_compat::par, local_hmatrix.get_UPLO(), local_hmatrix);
                            }
                            end = std::chrono::steady_clock::now();

                            std::chrono::duration<double> duration_facto = end - start;

                            // Timer for solve
                            start = std::chrono::steady_clock::now();
                            for (int i = 0; i < number_of_solves; i++) { // in place solve so Y_dense is variable
                                if (symmetry_type == 'N') {
                                    lu_solve(trans, local_hmatrix, Y_dense);
                                } else if (symmetry_type == 'S') {
                                    cholesky_solve(local_hmatrix.get_UPLO(), local_hmatrix, Y_dense);
                                }
                            }
                            end = std::chrono::steady_clock::now();

                            std::chrono::duration<double> duration_solve = end - start;

                            // data saving
                            savefile << epsilon << "," << dim << "," << n_threads << "," << policy_type << "," << id_rep << "," << compression_ratio << "," << space_saving << "," << duration_facto.count() << "," << duration_solve.count() << "," << clustering_type << "," << low_rank_generator_type << "," << symmetry_type << "," << generator_type << "," << hardware_type << "," << version << "\n";

                            List_factorization_duration[id_rep] = duration_facto.count();
                            List_solving_duration[id_rep]       = duration_solve.count();
                            list_compression_ratio[id_rep]      = compression_ratio;
                            list_space_saving[id_rep]           = space_saving;
                        } else {
                            std::cerr << "Unknown policy_type: " << policy_type << std::endl;
                            exit(1);
                        }
                    }
                    // mean and stddev saving
                    double mean_facto, std_dev_facto, mean_solve, std_dev_solve, mean_comp_ratio, std_dev_comp_ratio, mean_space_saving, std_dev_space_saving;

                    compute_standard_deviation(List_factorization_duration, number_of_repetitions, mean_facto, std_dev_facto);
                    compute_standard_deviation(List_solving_duration, number_of_repetitions, mean_solve, std_dev_solve);
                    compute_standard_deviation(list_compression_ratio, number_of_repetitions, mean_comp_ratio, std_dev_comp_ratio);
                    compute_standard_deviation(list_space_saving, number_of_repetitions, mean_space_saving, std_dev_space_saving);

                    savefile << epsilon << "," << dim << "," << n_threads << "," << policy_type << "," << "mean" << "," << mean_comp_ratio << "," << mean_space_saving << "," << mean_facto << "," << mean_solve << "," << clustering_type << "," << low_rank_generator_type << "," << symmetry_type << "," << generator_type << "," << hardware_type << "," << version << "\n";
                    savefile << epsilon << "," << dim << "," << n_threads << "," << policy_type << "," << "std_dev" << "," << std_dev_comp_ratio << "," << std_dev_space_saving << "," << std_dev_facto << "," << std_dev_solve << "," << clustering_type << "," << low_rank_generator_type << "," << symmetry_type << "," << generator_type << "," << hardware_type << "," << version << "\n";

                    is_ratio_done = true;
                }
            }
            id_pbl_size++;
        }
    }
    savefile.close();
}
#endif
