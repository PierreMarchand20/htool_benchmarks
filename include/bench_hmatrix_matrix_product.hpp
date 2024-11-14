#ifndef HTOOL_BENCHMARKS_HMATRIX_MATRIX_PRODUCT_HPP
#define HTOOL_BENCHMARKS_HMATRIX_MATRIX_PRODUCT_HPP

#include "NEW_add_hmatrix_vector_product.hpp" // où la fonction que l'on teste se situe

#include <htool/clustering/clustering.hpp>
#include <htool/hmatrix/hmatrix.hpp>
#include <htool/hmatrix/hmatrix_output.hpp>
#include <htool/hmatrix/linalg/interface.hpp>
#include <htool/hmatrix/lrmat/SVD.hpp>
#include <htool/hmatrix/lrmat/linalg/interface.hpp>
#include <htool/hmatrix/tree_builder/tree_builder.hpp>
#include <htool/testing/generate_test_case.hpp>
#include <htool/testing/generator_input.hpp>
#include <htool/testing/generator_test.hpp>
#include <htool/testing/geometry.hpp>
#include <htool/testing/partition.hpp>

#include "utils.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

/* Fixture */
class FT_LinearAlgebra {
  public:
    std::unique_ptr<htool::TestCaseProduct<double, htool::GeneratorTestDoubleSymmetric>> test_case;
    std::vector<double> B_vec, C_vec;
    double alpha, beta;
    std::unique_ptr<htool::HMatrix<double, htool::underlying_type<double>>> root_hmatrix;

    void SetUp(int n1, int n2, htool::underlying_type<double> epsilon, double eta, char transa) {
        test_case = std::make_unique<htool::TestCaseProduct<double, htool::GeneratorTestDoubleSymmetric>>(transa, 'N', n1, n2, 1, 1, 2, -1);

        const htool::Cluster<htool::underlying_type<double>> *root_cluster_A_output, *root_cluster_A_input, *root_cluster_B_output, *root_cluster_B_input, *root_cluster_C_output, *root_cluster_C_input;

        root_cluster_A_output = test_case->root_cluster_A_output;
        root_cluster_A_input  = test_case->root_cluster_A_input;
        root_cluster_B_output = test_case->root_cluster_B_output;
        // root_cluster_B_input  = &test_case->root_cluster_B_input->get_cluster_on_partition(0);
        root_cluster_C_output = test_case->root_cluster_C_output;
        // root_cluster_C_input  = &test_case->root_cluster_C_input->get_cluster_on_partition(0);

        // HMatrixTreeBuilder
        htool::HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*root_cluster_A_output, *root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*test_case->operator_A));

        // Dense matrix
        int no_B = root_cluster_B_output->get_size();
        int no_C = root_cluster_C_output->get_size();

        // Random Input matrix
        B_vec.resize(no_B);
        C_vec.resize(no_C);
        htool::generate_random_vector(B_vec);
        htool::generate_random_vector(C_vec);

        htool::generate_random_scalar(alpha);
        htool::generate_random_scalar(beta);
    }
};
#endif
