/// This program accepts an optional command-line argument to specify the symmetry type
/// of the H-matrix. The default symmetry type is 'N' (non-symmetric). The program measures
/// and outputs the time taken to perform the benchmark.

#include "htool_benchmark/bench_hmatrix_factorization.hpp"
#include "htool_benchmark/generator_types.hpp"
#include "htool_benchmark/hmatrix_fixture.hpp"

using namespace htool_benchmark;

int main(int argc, char *argv[]) {
    // Check if the number of arguments is valid
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " <symmetry_type>[default=N] = {N|S}" << std::endl;
        return 1;
    }

    // Determine the symmetry type; default to 'N' if not provided
    const char symmetry_type = (argc == 2) ? std::toupper(argv[1][0]) : 'N';

    // Validate the symmetry type
    if (symmetry_type != 'N' && symmetry_type != 'S') {
        std::cerr << "Error: invalid symmetry type. Must be either 'N' or 'S'." << std::endl;
        return 1;
    }

    // Measure the time taken to perform the benchmark
    auto start = std::chrono::high_resolution_clock::now();
    bench_hmatrix_factorization<FixtureHMatrix<FixtureGenerator<LaplaceLikeGenerator>>>(symmetry_type);
    auto end                           = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    // Output the time taken for the benchmark
    std::cout << "Time taken: " << diff.count() << " seconds" << std::endl;
    std::cout << "================================================" << std::endl
              << std::endl;

    return 0;
}
