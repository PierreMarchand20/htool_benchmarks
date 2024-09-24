#include "bench_hmatrix_matrix_product.hpp"
#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
    const int number_of_repetitions = 9;
    const int min_dim_pbl           = 1 << 16; // >= 1 << 16 on laptop else variance too high
    const int max_dim_pbl           = 1 << 18; // <= (1 << 15) for dense or (1 << 18) otherwise on laptop else "Abandon (core dumped)"
    const int dim_pbl_step          = 2;

    std::ofstream savefile;
    savefile.open("bench_hmatrix_matrix_product_vs_pbl_size.csv");
    savefile << "epsilon, dim_pbl, algo_type, id_rep, compression_ratio, space_saving, time (s) | mean time (s) | standard_deviation \n";

    for (double epsilon : {1e-10, 1e-8, 1e-6}) {
        for (int dim_pbl = min_dim_pbl; dim_pbl <= max_dim_pbl; dim_pbl *= dim_pbl_step) {
            // Setup
            FT_LinearAlgebra Fixture;
            double eta  = 10;
            char transa = 'N';
            Fixture.SetUp(dim_pbl, dim_pbl, epsilon, eta, transa);
            double List_duration[number_of_repetitions] = {0};

            for (string algo_type : {"Classic", "TaskBased"}) { // max_dim_pbl <= (1 << 15) for dense on laptop else "Abandon (core dumped)"
                for (int id_rep = 0; id_rep < number_of_repetitions; id_rep++) {
                    std::chrono::steady_clock::time_point start, end;
                    double compression_ratio = 0;
                    double space_saving      = 0;
                    if (algo_type == "Classic") {
                        // Timer
                        start = std::chrono::steady_clock::now();
                        openmp_internal_add_hmatrix_vector_product(transa, Fixture.alpha, *Fixture.root_hmatrix, Fixture.B_vec.data(), Fixture.beta, Fixture.C_vec.data());
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
                        NEW_openmp_add_hmatrix_vector_product(transa, Fixture.alpha, *Fixture.root_hmatrix, Fixture.B_vec.data(), Fixture.beta, Fixture.C_vec.data());
                        end                                    = std::chrono::steady_clock::now();
                        std::chrono::duration<double> duration = end - start;

                        // Compression ratio and space saving
                        auto hmatrix_information = get_hmatrix_information(*Fixture.root_hmatrix);
                        compression_ratio        = std::stod(hmatrix_information["Compression_ratio"]);
                        space_saving             = std::stod(hmatrix_information["Space_saving"]);
                        // std::cout << hmatrix_information["Number_of_threads"] << std::endl;

                        // data saving
                        savefile << epsilon << ", " << dim_pbl << ", " << algo_type << ", " << id_rep << ", " << compression_ratio << ", " << space_saving << ", " << duration.count() << "\n";
                        List_duration[id_rep] = duration.count();
                    } else if (algo_type == "Dense") { // max_dim_pbl <= (1 << 15) for dense on laptop else "Abandon (core dumped)"
                        std::unique_ptr<Matrix<double>> HA_dense;
                        HA_dense = std::make_unique<Matrix<double>>(Fixture.root_hmatrix->get_target_cluster().get_size(), Fixture.root_hmatrix->get_source_cluster().get_size());
                        copy_to_dense(*Fixture.root_hmatrix, HA_dense->data());

                        // Timer
                        start = std::chrono::steady_clock::now();
                        add_matrix_vector_product(transa, Fixture.alpha, *HA_dense, Fixture.B_vec.data(), Fixture.beta, Fixture.C_vec.data());
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