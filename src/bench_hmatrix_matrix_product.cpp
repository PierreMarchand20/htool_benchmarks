#include "htool_benchmark/bench_hmatrix_matrix_product.hpp"
#include "htool_benchmark/cli.hpp"
#include "htool_benchmark/generator_fixture.hpp"
#include "htool_benchmark/generator_types.hpp"
#include "htool_benchmark/hmatrix_fixture.hpp"

using namespace htool_benchmark;

int main(int argc, char *argv[]) {

    bool is_error = false;

    is_error = is_error || check_input_size(argc, argv);

    // Convert the first argument to an integer test case id
    std::string test_case_type;
    is_error = is_error || check_test_case_type(argc, argv, test_case_type);

    // Determine the symmetry type; default to 'N' if not provided
    char symmetry_type;
    is_error = is_error || check_symmetry_type(argc, argv, symmetry_type);

    //
    std::string generator_type;
    is_error = is_error || check_generator_type(argc, argv, generator_type);

    //
    std::string clustering_type;
    is_error = is_error || check_clustering_type(argc, argv, clustering_type);

    //
    std::string low_rank_generator_type;
    is_error = is_error || check_low_rank_generator_type(argc, argv, low_rank_generator_type);

    //
    std::string hardware_type;
    is_error = is_error || check_hardware_type(argc, argv, hardware_type);

    //
    std::string version;
    is_error = is_error || check_version(argc, argv, version);

    if (is_error)
        return 1;

    // Measure the time taken to perform the benchmark
    auto start = std::chrono::high_resolution_clock::now();
    if (generator_type == "Laplace")
        bench_hmatrix_matrix_product<FixtureHMatrix, FixtureGenerator, OptimizedLaplaceLikeGenerator>(test_case_type, symmetry_type, generator_type, clustering_type, low_rank_generator_type, hardware_type, version);
    else if (generator_type == "Helmholtz")
        bench_hmatrix_matrix_product<FixtureHMatrix, FixtureGenerator, OptimizedHelmholtzLikeGenerator>(test_case_type, symmetry_type, generator_type, clustering_type, low_rank_generator_type, hardware_type, version);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;

    // Output the time taken
    std::cout << "Time taken: " << diff.count() << " seconds" << std::endl;
    std::cout << "================================================" << std::endl
              << std::endl;

    return 0;
}
