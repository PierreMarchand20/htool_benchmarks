/// This program benchmarks the H-matrix build process using the
/// specified test case and symmetry type. It takes two command line
/// arguments: the first is the id of the test case to run, and the
/// second is the symmetry type of the test case. If the second
/// argument is not provided, it defaults to 'N' (no symmetry).

#include "htool_benchmark/bench_hmatrix_build.hpp"
#include "htool_benchmark/generator_fixture.hpp"
#include "htool_benchmark/generator_types.hpp"
#include "utils.hpp"
#include <string>

using namespace htool_benchmark;

int main(int argc, char *argv[]) {
    // Check for valid number of arguments
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <test_case> = {0|1|2} <symmetry_type>[default=N] = {N|S}" << std::endl;
        return 1;
    }

    // Convert the first argument to an integer test case id
    const int test_case_id = std::stoi(argv[1]);
    // Get the name of the test case based on the provided id
    const std::string test_case_name = get_test_case_name(test_case_id);

    // Determine the symmetry type; default to 'N' if not provided
    const char symmetry_type = (argc == 3) ? std::toupper(argv[2][0]) : 'N';

    // Validate symmetry type
    if (symmetry_type != 'N' && symmetry_type != 'S') {
        std::cerr << "Error: invalid symmetry type. Must be either 'N' or 'S'." << std::endl;
        return 1;
    }

    // Measure the time taken to perform the benchmark
    auto start = std::chrono::high_resolution_clock::now();
    bench_hmatrix_build<FixtureGenerator<LaplaceLikeGenerator>>(test_case_name, symmetry_type);
    auto end                           = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    // Output the time taken
    std::cout << "Time taken: " << diff.count() << " seconds" << std::endl;
    std::cout << "================================================" << std::endl
              << std::endl;

    return 0;
}
