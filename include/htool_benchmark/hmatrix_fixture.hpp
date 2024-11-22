#include "generator_fixture.hpp"

template <typename FixtureGenerator>
class FixtureHMatrix {

    std::unique_ptr<FixtureGenerator> fixture_generator;
    std::unique_ptr<htool::HMatrix<double, htool::underlying_type<double>>> root_hmatrix;

  public:
    void setup_benchmark(int N, htool::underlying_type<double> epsilon, double eta, char symmetry) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(N);
        htool::HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*fixture_generator->m_target_root_cluster, *fixture_generator->m_target_root_cluster, epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L', -1, -1, -1);
        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(fixture_generator->generator.get()));
    }

    void setup_benchmark(int M, int N, htool::underlying_type<double> epsilon, double eta) {
        fixture_generator = std::make_unique<FixtureGenerator>();
        fixture_generator->setup_benchmark(M, N);
        htool::HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*fixture_generator->m_target_root_cluster, *fixture_generator->m_source_root_cluster, epsilon, eta, 'N', 'N', -1, -1, -1);
        root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(fixture_generator->generator.get()));
    }

    // std::unique_ptr<htool::TestCaseProduct<double, htool::GeneratorTestDoubleSymmetric>> test_case;
    // std::vector<double> B_vec, C_vec;
    // double alpha, beta;

    // void SetUp(int n1, int n2, htool::underlying_type<double> epsilon, double eta, char transa) {
    //     test_case = std::make_unique<htool::TestCaseProduct<double, htool::GeneratorTestDoubleSymmetric>>(transa, 'N', n1, n2, 1, 1, 2, -1);

    //     const htool::Cluster<htool::underlying_type<double>> *root_cluster_A_output, *root_cluster_A_input, *root_cluster_B_output, *root_cluster_B_input, *root_cluster_C_output, *root_cluster_C_input;

    //     root_cluster_A_output = test_case->root_cluster_A_output;
    //     root_cluster_A_input  = test_case->root_cluster_A_input;
    //     root_cluster_B_output = test_case->root_cluster_B_output;
    //     // root_cluster_B_input  = &test_case->root_cluster_B_input->get_cluster_on_partition(0);
    //     root_cluster_C_output = test_case->root_cluster_C_output;
    //     // root_cluster_C_input  = &test_case->root_cluster_C_input->get_cluster_on_partition(0);

    //     // HMatrixTreeBuilder
    //     htool::HMatrixTreeBuilder<double, htool::underlying_type<double>> hmatrix_tree_builder(*root_cluster_A_output, *root_cluster_A_input, epsilon, eta, 'N', 'N', -1, -1, -1);
    //     root_hmatrix = std::make_unique<htool::HMatrix<double, htool::underlying_type<double>>>(hmatrix_tree_builder.build(*test_case->operator_A));

    //     // Dense matrix
    //     int no_B = root_cluster_B_output->get_size();
    //     int no_C = root_cluster_C_output->get_size();

    //     // Random Input matrix
    //     B_vec.resize(no_B);
    //     C_vec.resize(no_C);
    //     htool::generate_random_vector(B_vec);
    //     htool::generate_random_vector(C_vec);

    //     htool::generate_random_scalar(alpha);
    //     htool::generate_random_scalar(beta);
    // }
};