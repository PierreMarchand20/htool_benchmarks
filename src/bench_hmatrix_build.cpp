#include "htool_benchmark/bench_hmatrix_build.hpp"
#include "htool_benchmark/cli.hpp"
#include "htool_benchmark/generator_fixture.hpp"
#include "htool_benchmark/generator_types.hpp"
#include <string>

using namespace htool_benchmark;

int main(int argc, char *argv[]) {

    bool is_error = false;

    is_error = is_error || check_input_size(argc, argv);

    //
    std::string test_case_type;
    is_error = is_error || check_test_case_type(argc, argv, test_case_type);

    //
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
    std::string policy_type;
    is_error = is_error || check_policy_type(argc, argv, policy_type);

    if (is_error)
        return 1;

    // Measure the time taken to perform the benchmark
    auto start = std::chrono::high_resolution_clock::now();
    if (generator_type == "Laplace")
        bench_hmatrix_build<FixtureGenerator, OptimizedLaplaceLikeGenerator>(test_case_type, symmetry_type, generator_type, clustering_type, low_rank_generator_type, policy_type);
    else if (generator_type == "Helmholtz")
        bench_hmatrix_build<FixtureGenerator, OptimizedHelmholtzLikeGenerator>(test_case_type, symmetry_type, generator_type, clustering_type, low_rank_generator_type, policy_type);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;

    // Output the time taken
    std::cout << "Time taken: " << diff.count() << " seconds" << std::endl;
    std::cout << "================================================" << std::endl
              << std::endl;

    return 0;
}
