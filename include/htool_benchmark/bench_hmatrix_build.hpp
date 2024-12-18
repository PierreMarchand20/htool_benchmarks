#ifndef HTOOL_BENCHMARK_HMATRIX_BUILD_HPP
#define HTOOL_BENCHMARK_HMATRIX_BUILD_HPP

#include "NEW_tree_builder.hpp"
#include "generator_fixture.hpp"
#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
#include <iostream>

using namespace htool;

namespace htool_benchmark {

template <typename FixtureGenerator>
void bench_hmatrix_build(std::string test_case_type, char symmetry_type) {

    // declare variables
    std::vector<std::string> List_algo_type;
    std::vector<double> List_epsilon;
    std::vector<int> List_pbl_size;
    std::vector<int> List_thread;
    int id_pbl_size(0), id_thread(0);
    bool is_ratio_done(false);

    // custom parameters
    const int number_of_repetitions = 2;
    List_algo_type                  = {"Classic", "TaskBased"};
    List_epsilon                    = {1e-10, 1e-8, 1e-6, 1e-4};

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
        List_thread   = {1, 2, 4, 8, 16};
    }

    // header csv file
    std::ofstream savefile;
    savefile.open("bench_hmatrix_build_vs_" + test_case_type + ".csv");
    savefile << "epsilon, dim_pbl, number_of_threads, algo_type, id_rep, compression_ratio, space_saving, time (s) \n";

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
            FixtureGenerator fixture;
            double eta = 10;
            const htool::Cluster<double> *target_cluster, *source_cluster;
            if (symmetry_type != 'N') {
                fixture.setup_benchmark(dim_pbl);
                target_cluster = fixture.m_target_root_cluster.get();
                source_cluster = target_cluster;
            } else {
                fixture.setup_benchmark(dim_pbl, dim_pbl);
                target_cluster = fixture.m_target_root_cluster.get();
                source_cluster = fixture.m_source_root_cluster.get();
            }
            double list_build_duration[number_of_repetitions]    = {0};
            double list_compression_ratio[number_of_repetitions] = {0};
            double list_space_saving[number_of_repetitions]      = {0};

            for (std::string algo_type : List_algo_type) {
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
                        double compression_ratio = 0;
                        double space_saving      = 0;
                        if (algo_type == "Classic") {
                            // Hmatrix
                            using HMatrixTreeBuilderType = htool::HMatrixTreeBuilder<double, htool::underlying_type<double>>;
                            HMatrixTreeBuilderType hmatrix_tree_builder(*target_cluster, *source_cluster, epsilon, eta, symmetry_type, symmetry_type == 'N' ? 'N' : 'L', -1, -1, -1);

                            // Timer
                            start             = std::chrono::steady_clock::now();
                            auto root_hmatrix = hmatrix_tree_builder.build(*fixture.generator); // *generator
                            end               = std::chrono::steady_clock::now();

                            std::chrono::duration<double> classic_build_duration = end - start;

                            // Compression ratio and space saving
                            auto hmatrix_information = get_hmatrix_information(root_hmatrix);
                            compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                            space_saving             = std::stod(hmatrix_information["Space_saving"]);

                            // data saving
                            savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << classic_build_duration.count() << "\n";
                            list_build_duration[id_rep]    = classic_build_duration.count();
                            list_compression_ratio[id_rep] = compression_ratio;
                            list_space_saving[id_rep]      = space_saving;

                        } else if (algo_type == "TaskBased") {
                            using HMatrixTreeBuilderType = htool::HMatrixTaskBasedTreeBuilder<double, htool::underlying_type<double>>;
                            HMatrixTreeBuilderType hmatrix_tree_builder(*target_cluster, *source_cluster, epsilon, eta, symmetry_type, symmetry_type == 'N' ? 'N' : 'L', -1, -1, -1);

                            // Timer
                            start             = std::chrono::steady_clock::now();
                            auto root_hmatrix = hmatrix_tree_builder.build(*fixture.generator); // *generator
                            end               = std::chrono::steady_clock::now();

                            std::chrono::duration<double> task_based_build_duration = end - start;

                            // Compression ratio and space saving
                            auto hmatrix_information = get_hmatrix_information(root_hmatrix);
                            compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                            space_saving             = std::stod(hmatrix_information["Space_saving"]);

                            // data saving
                            savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << task_based_build_duration.count() << "\n";
                            list_build_duration[id_rep]    = task_based_build_duration.count();
                            list_compression_ratio[id_rep] = compression_ratio;
                            list_space_saving[id_rep]      = space_saving;
                        }
                    }
                    // mean and stddev saving
                    double mean_build, std_dev_build, mean_comp_ratio, std_dev_comp_ratio, mean_space_saving, std_dev_space_saving;

                    compute_standard_deviation(list_build_duration, number_of_repetitions, mean_build, std_dev_build);
                    compute_standard_deviation(list_compression_ratio, number_of_repetitions, mean_comp_ratio, std_dev_comp_ratio);
                    compute_standard_deviation(list_space_saving, number_of_repetitions, mean_space_saving, std_dev_space_saving);

                    savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << "mean" << ", " << mean_comp_ratio << ", " << mean_space_saving << ", " << mean_build << "\n";
                    savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << "stddev" << ", " << std_dev_comp_ratio << ", " << std_dev_space_saving << ", " << std_dev_build << "\n";

                    is_ratio_done = true;
                }
            }
            id_pbl_size++;
        }
    }
    savefile.close();
}
} // namespace htool_benchmark
#endif