#ifndef HTOOL_BENCHMARK_HMATRIX_BUILD_HPP
#define HTOOL_BENCHMARK_HMATRIX_BUILD_HPP

#include "htool/hmatrix/execution_policies.hpp"
#include "htool_benchmark/cli.hpp"
#include "utils.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <htool/clustering/tree_builder/tree_builder.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
#include <iostream>
#include <string>

using namespace htool;

namespace htool_benchmark {
/**
 * @brief Benchmark of the H-matrix build process
 * @tparam FixtureGenerator The fixture used for setting up the generator.
 * @param test_case_type Type of the test case: "pbl_size","thread" or "ratio"
 * @param symmetry_type Symmetry type of the H-matrix: 'N' for non-symmetric, 'S' for symmetric, 'H' for hermitian
 *
 * This function benchmarks the H-matrix build process using the specified test case and symmetry type.
 *
 * The function measures the time taken to build the H-matrix using the Classic and TaskBased algorithms.
 * It also measures the compression ratio and space saving of the H-matrix.
 * The results are saved in a CSV file that can be read by plot_bench_vs_*.py.
 */
template <template <typename> class FixtureGenerator, typename GeneratorType>
void bench_hmatrix_build(std::string test_case_type, char symmetry_type, std::string generator_type, std::string clustering_type, std::string low_rank_generator_type, std::string policy_type, std::string hardware_type, std::string version) {
    using CoefficientPrecision = typename FixtureGenerator<GeneratorType>::CoefficientPrecision;

    // declare variables
    std::vector<double> List_epsilon;
    std::vector<int> List_pbl_size;
    std::vector<int> List_thread;
    int id_pbl_size(0), id_thread(0);
    bool is_ratio_done(false);

    // custom parameters
    const int number_of_repetitions = 1;
    List_epsilon                    = {1e-10, 1e-7, 1e-4};
    double eta                      = 10;

    if (test_case_type == "pbl_size") { // 1<<19 vs 1 thread OK sur Cholesky, 1<<20 vs 1 thread out of memory
        // List_pbl_size = {1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19};
        List_pbl_size = {1 << 8, 1 << 9};
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
    std::string filename = "bench_hmatrix_build_vs_" + test_case_type + ".csv";
    std::ofstream savefile;
    bool file_already_exists = std::filesystem::exists(filename);
    savefile.open(filename, std::ios::app);
    if (!file_already_exists) {
        savefile << "epsilon,size,number_of_threads,policy_type,id_rep,compression_ratio,space_saving,time, clustering_type,low_rank_generator_type,symmetry_type,generator_type,hardware_type,version\n";
    }

    // cout parameters
    std::cout << "++++++++++++++++++ Test case: ++++++++++++++++++" << std::endl;
    std::cout << "Test_case: " << test_case_type << std::endl;
    std::cout << "List_epsilon: " << List_epsilon << std::endl;
    std::cout << "List_pbl_size: " << List_pbl_size << std::endl;
    std::cout << "List_thread: " << List_thread << std::endl;
    std::cout << "Number_of_repetitions: " << number_of_repetitions << std::endl;
    std::cout << "Symmetry_type: " << symmetry_type << std::endl;
    std::cout << "Generator_type: " << generator_type << std::endl;
    std::cout << "Clustering_type: " << clustering_type << std::endl;
    std::cout << "Low_rank_generator_type: " << low_rank_generator_type << std::endl;
    std::cout << "Policy_type: " << policy_type << std::endl;
    std::cout << "Hardware_type: " << hardware_type << std::endl;
    std::cout << "Version: " << version << std::endl;
    std::cout << "Eta: " << eta << std::endl;
    std::cout << std::endl;

    // Partitioning strategy
    auto partitioning_strategy = process_clustering_type<double>(clustering_type);
    auto cluster_tree_builder  = std::make_shared<htool::ClusterTreeBuilder<double>>();
    cluster_tree_builder->set_partitioning_strategy(partitioning_strategy);

    // computation
    for (double epsilon : List_epsilon) {
        id_pbl_size = 0;

        for (int size : List_pbl_size) {
            // Setup
            FixtureGenerator<GeneratorType> fixture(cluster_tree_builder);
            const htool::Cluster<double> *target_cluster, *source_cluster;

            fixture.setup_benchmark(size);
            target_cluster = fixture.m_target_root_cluster.get();
            source_cluster = fixture.m_source_root_cluster.get();

            HMatrixTreeBuilder<CoefficientPrecision> hmatrix_tree_builder(epsilon, eta, symmetry_type, symmetry_type == 'N' ? 'N' : 'L');
            std::shared_ptr<htool::VirtualInternalLowRankGenerator<CoefficientPrecision>> compression_strategy;
            if constexpr (GeneratorType::require_permuted_input) {
                compression_strategy = process_low_rank_generator_type<CoefficientPrecision>(low_rank_generator_type, *fixture.generator);
            } else {
                compression_strategy = process_low_rank_generator_type<CoefficientPrecision>(low_rank_generator_type, *fixture.generator, target_cluster->get_permutation(), source_cluster->get_permutation());
            }
            hmatrix_tree_builder.set_low_rank_generator(compression_strategy);

            double list_build_duration[number_of_repetitions]    = {0};
            double list_compression_ratio[number_of_repetitions] = {0};
            double list_space_saving[number_of_repetitions]      = {0};

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
                    double compression_ratio;
                    double space_saving;
                    std::chrono::duration<double> duration;

                    if (policy_type == "seq") {
                        // Timer
                        start             = std::chrono::steady_clock::now();
                        auto root_hmatrix = hmatrix_tree_builder.build(exec_compat::seq, *fixture.generator, *target_cluster, *source_cluster);
                        end               = std::chrono::steady_clock::now();

                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(root_hmatrix);
                        compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                        space_saving             = std::stod(hmatrix_information["Space_saving"]);
                    } else if (policy_type == "par") {
                        // Timer
                        start             = std::chrono::steady_clock::now();
                        auto root_hmatrix = hmatrix_tree_builder.build(exec_compat::par, *fixture.generator, *target_cluster, *source_cluster);
                        end               = std::chrono::steady_clock::now();

                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(root_hmatrix);
                        compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                        space_saving             = std::stod(hmatrix_information["Space_saving"]);

                    } else if (policy_type == "omp_task") {
                        omp_task_policy<CoefficientPrecision, double> omp_task;
                        // Timer
                        start             = std::chrono::steady_clock::now();
                        auto root_hmatrix = hmatrix_tree_builder.build(omp_task, *fixture.generator, *target_cluster, *source_cluster);
                        end               = std::chrono::steady_clock::now();

                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(root_hmatrix);
                        compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                        space_saving             = std::stod(hmatrix_information["Space_saving"]);
                    } else {
                        std::cerr << "Unknown policy_type: " << policy_type << std::endl;
                        exit(1);
                    }

                    duration = end - start;

                    // data saving
                    savefile << epsilon << "," << size << "," << n_threads << "," << policy_type << "," << id_rep << "," << compression_ratio << "," << space_saving << "," << duration.count() << "," << clustering_type << "," << low_rank_generator_type << "," << symmetry_type << "," << generator_type << "," << hardware_type << "," << version << "\n";
                    list_build_duration[id_rep]    = duration.count();
                    list_compression_ratio[id_rep] = compression_ratio;
                    list_space_saving[id_rep]      = space_saving;
                }
                // mean and stddev saving
                double mean_build, std_dev_build, mean_comp_ratio, std_dev_comp_ratio, mean_space_saving, std_dev_space_saving;

                compute_standard_deviation(list_build_duration, number_of_repetitions, mean_build, std_dev_build);
                compute_standard_deviation(list_compression_ratio, number_of_repetitions, mean_comp_ratio, std_dev_comp_ratio);
                compute_standard_deviation(list_space_saving, number_of_repetitions, mean_space_saving, std_dev_space_saving);

                savefile << epsilon << "," << size << "," << n_threads << "," << policy_type << "," << "mean" << "," << mean_comp_ratio << "," << mean_space_saving << "," << mean_build << "," << clustering_type << "," << low_rank_generator_type << "," << symmetry_type << "," << generator_type << "," << hardware_type << "," << version << "\n";
                savefile << epsilon << "," << size << "," << n_threads << "," << policy_type << "," << "stddev" << "," << std_dev_comp_ratio << "," << std_dev_space_saving << "," << std_dev_build << "," << clustering_type << "," << low_rank_generator_type << "," << symmetry_type << "," << generator_type << "," << hardware_type << "," << version << "\n";

                is_ratio_done = true;
            }

            id_pbl_size++;
        }
    }
    savefile.close();
}
} // namespace htool_benchmark
#endif
