#include "benchmark_hmatrix_build.hpp"

using namespace std;
using namespace htool;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    bool is_error            = false;
    const int number_of_rows = 200;
    // const int number_of_rows_increased    = 400;
    const int number_of_columns = 200;
    // const int number_of_columns_increased = 400;

    // for (auto use_local_cluster : {true, false}) {
    //     for (auto epsilon : {1e-14, 1e-6}) {
    //         for (auto use_dense_block_generator : {true, false}) {
    //             // Square matrix
    //             is_error = is_error || test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, use_local_cluster, 'N', 'N', epsilon, use_dense_block_generator);
    //             is_error = is_error || test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, use_local_cluster, 'S', 'U', epsilon, use_dense_block_generator);
    //             is_error = is_error || test_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, use_local_cluster, 'S', 'L', epsilon, use_dense_block_generator);

    //             // Rectangle matrix
    //             is_error = is_error || test_hmatrix_build<double, GeneratorTestDouble>(number_of_rows_increased, number_of_columns, use_local_cluster, 'N', 'N', epsilon, use_dense_block_generator);
    //             is_error = is_error || test_hmatrix_build<double, GeneratorTestDouble>(number_of_rows, number_of_columns_increased, use_local_cluster, 'N', 'N', epsilon, use_dense_block_generator);
    //         }
    //     }
    // }
    is_error = is_error || benchmark_hmatrix_build<double, GeneratorTestDoubleSymmetric>(number_of_rows, number_of_columns, true, 'N', 'N', 1e-14);

    MPI_Finalize();

    if (is_error) {
        return 1;
    }
    return 0;
}
