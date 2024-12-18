#include "htool_benchmark/bench_hmatrix_factorization.hpp"
#include "htool_benchmark/generator_types.hpp"
#include "htool_benchmark/hmatrix_fixture.hpp"

using namespace htool_benchmark;

int main(int argc, char *argv[]) {
    // bench_hmatrix_factorization<double, GeneratorTestDouble>();
    bench_hmatrix_factorization<FixtureHMatrix<FixtureGenerator<LaplaceLikeGenerator>>>('N');
    return 0;
}
