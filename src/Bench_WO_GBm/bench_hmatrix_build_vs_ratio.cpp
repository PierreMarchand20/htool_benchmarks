#include "bench_hmatrix_build.hpp"
#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std;
using namespace htool;

// #include <vector>
// #include <cmath>

std::vector<int> geometric_progression(int a, int r, int n) {
    // Computes the sequence of n terms of a geometric progression.
    //
    // Parameters:
    // a: The first term of the progression.
    // r: The common ratio of the progression.
    // n: The number of terms to compute.
    //
    // Returns:
    // A vector containing the sequence of n terms of the progression.
    std::vector<int> sequence(n);
    sequence[0] = a;
    for (int i = 1; i < n; i++) {
        sequence[i] = sequence[i - 1] * r;
    }
    return sequence;
}

int main(int argc, char *argv[]) {
    const int number_of_repetitions = 2;
    const int number_of_cases       = 5;
    std::vector<int> List_pbl_size(number_of_cases);
    std::vector<int> List_thread(number_of_cases);
    List_pbl_size = geometric_progression(1 << 14, 2, number_of_cases);
    List_thread   = geometric_progression(1, 2, number_of_cases);

    std::ofstream savefile;
    savefile.open("bench_hmatrix_build_vs_ratio.csv");
    savefile << "epsilon, dim_pbl, number_of_threads, algo_type, id_rep, compression_ratio, space_saving, time (s) | mean time (s) | standard_deviation \n";

    for (double epsilon : {1e-10, 1e-8, 1e-6}) {
        int id_thread = 0;
        for (int dim_pbl : List_pbl_size) {
            // Setup
            FT_Generator fixture;
            double eta = 10;
            fixture.SetUp(dim_pbl, dim_pbl, 'N', 'N', epsilon, eta);
            double list_build_duration[number_of_repetitions] = {0};

            int n_threads = List_thread[id_thread];
            id_thread++;

            for (string algo_type : {"Classic", "TaskBased"}) {
                for (int id_rep = 0; id_rep < number_of_repetitions; id_rep++) {
                    std::chrono::steady_clock::time_point start, end;
                    omp_set_num_threads(n_threads);

                    if (algo_type == "Classic" || algo_type == "TaskBased") { // Hmatrix
                        double compression_ratio = 0;
                        double space_saving      = 0;
                        if (algo_type == "Classic") {
                            // Hmatrix
                            using HMatrixTreeBuilderType = HMatrixTreeBuilder<double, htool::underlying_type<double>>;
                            std::unique_ptr<HMatrixTreeBuilderType> hmatrix_tree_builder;
                            hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilderType>(fixture.m_target_root_cluster->get_cluster_on_partition(0), fixture.m_source_root_cluster->get_cluster_on_partition(0), epsilon, eta, 'N', 'N', -1, -1, -1);

                            // Timer
                            start                                                = std::chrono::steady_clock::now();
                            auto root_hmatrix                                    = hmatrix_tree_builder->build(*fixture.generator); // *generator
                            end                                                  = std::chrono::steady_clock::now();
                            std::chrono::duration<double> classic_build_duration = end - start;

                            // Compression ratio and space saving
                            auto hmatrix_information = get_hmatrix_information(root_hmatrix);
                            compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                            space_saving             = std::stod(hmatrix_information["Space_saving"]);

                            // data saving
                            savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << classic_build_duration.count() << "\n";
                            list_build_duration[id_rep] = classic_build_duration.count();

                        } else if (algo_type == "TaskBased") {
                            using HMatrixTreeBuilderType = HMatrixTaskBasedTreeBuilder<double, htool::underlying_type<double>>;
                            std::unique_ptr<HMatrixTreeBuilderType> hmatrix_tree_builder;
                            hmatrix_tree_builder = std::make_unique<HMatrixTreeBuilderType>(fixture.m_target_root_cluster->get_cluster_on_partition(0), fixture.m_source_root_cluster->get_cluster_on_partition(0), epsilon, eta, 'N', 'N', -1, -1, -1);

                            // Timer
                            start                                                   = std::chrono::steady_clock::now();
                            auto root_hmatrix                                       = hmatrix_tree_builder->build(*fixture.generator); // *generator
                            end                                                     = std::chrono::steady_clock::now();
                            std::chrono::duration<double> task_based_build_duration = end - start;

                            // Compression ratio and space saving
                            auto hmatrix_information = get_hmatrix_information(root_hmatrix);
                            compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                            space_saving             = std::stod(hmatrix_information["Space_saving"]);

                            // data saving
                            savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << task_based_build_duration.count() << "\n";
                            list_build_duration[id_rep] = task_based_build_duration.count();
                        }
                        // } else { // Dense. Max dimension pbl = 1 << 15
                        //     std::vector<double> dense_data(dim_pbl * dim_pbl);

                        //     // Timer
                        //     start = std::chrono::steady_clock::now();
                        //     (*fixture.generator).copy_submatrix(dim_pbl, dim_pbl, 0, 0, dense_data.data());
                        //     end                                                = std::chrono::steady_clock::now();
                        //     std::chrono::duration<double> dense_build_duration = end - start;

                        //     // data saving
                        //     savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << id_rep << ", " << 0 << ", " << 0 << ", " << dense_build_duration.count() << "\n";
                        //     list_build_duration[id_rep] = dense_build_duration.count();
                    }
                }
                // mean and stddev saving
                double mean, std_dev;
                compute_standard_deviation(list_build_duration, number_of_repetitions, mean, std_dev);
                savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << "mean" << ", " << "N.A" << ", " << "N.A" << ", " << mean << "\n";
                savefile << epsilon << ", " << dim_pbl << ", " << n_threads << ", " << algo_type << ", " << "stddev" << ", " << "N.A" << ", " << "N.A" << ", " << std_dev << "\n";
            }
        }
    }
    savefile.close();
    return 0;
}
