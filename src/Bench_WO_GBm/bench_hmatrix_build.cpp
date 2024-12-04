#include "htool_benchmark/bench_hmatrix_build.hpp"
#include "htool_benchmark/generator_fixture.hpp"
#include "htool_benchmark/generator_types.hpp"

using namespace htool_benchmark;

int main(int argc, char *argv[]) {
    // get test_case id
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <test_case> " << std::endl;
        return 1;
    }
    int test_case = std::stoi(argv[1]); // 0: pbl_size, 1: threads, 2: ratio pbl_size/thread

    std::vector<std::string> Dict_test_case = {"pbl_size", "thread", "ratio"};
    std::string case_name                   = Dict_test_case[test_case];

    bench_hmatrix_build<FixtureGenerator<LaplaceLikeGenerator>>(case_name, 'S');
    return 0;
}
