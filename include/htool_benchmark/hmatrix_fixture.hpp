#ifndef HTOOL_BENCHMARK_HMATRIX_FIXTURE_HPP
#define HTOOL_BENCHMARK_HMATRIX_FIXTURE_HPP

#include "cli.hpp"
#include "htool/clustering/cluster_node.hpp"
#include <htool/hmatrix/tree_builder/tree_builder.hpp>

using namespace std;
using namespace htool;
template <template <typename> class FixtureGenerator, typename GeneratorType>
class FixtureHMatrix {
    using CoefficientPrecision = typename FixtureGenerator<GeneratorType>::CoefficientPrecision;
    std::shared_ptr<FixtureGenerator<GeneratorType>> m_fixture_generator;

  public:
    std::unique_ptr<htool::HMatrix<CoefficientPrecision, htool::underlying_type<CoefficientPrecision>>> root_hmatrix;
    std::vector<HMatrix<CoefficientPrecision> *> L0;
    int max_node_L0 = 64;

    FixtureHMatrix(std::shared_ptr<FixtureGenerator<GeneratorType>> fixture_generator) : m_fixture_generator(fixture_generator) {}

    void setup_benchmark(int N, htool::underlying_type<CoefficientPrecision> epsilon, double eta, char symmetry, std::string low_rank_generator_type) {
        m_fixture_generator->setup_benchmark(N);
        const htool::Cluster<double> *target_cluster, *source_cluster;
        target_cluster = m_fixture_generator->m_target_root_cluster.get();
        source_cluster = m_fixture_generator->m_source_root_cluster.get();

        HMatrixTreeBuilder<CoefficientPrecision> hmatrix_tree_builder(epsilon, eta, symmetry, symmetry == 'N' ? 'N' : 'L');

        std::shared_ptr<htool::VirtualInternalLowRankGenerator<CoefficientPrecision>> compression_strategy;
        if constexpr (GeneratorType::require_permuted_input) {
            compression_strategy = process_low_rank_generator_type<CoefficientPrecision>(low_rank_generator_type, *m_fixture_generator->generator);
        } else {
            compression_strategy = process_low_rank_generator_type<CoefficientPrecision>(low_rank_generator_type, *m_fixture_generator->generator, target_cluster->get_permutation(), source_cluster->get_permutation());
        }
        hmatrix_tree_builder.set_low_rank_generator(compression_strategy);

        root_hmatrix = std::make_unique<htool::HMatrix<CoefficientPrecision, double>>(hmatrix_tree_builder.openmp_build(*(m_fixture_generator->generator), *(target_cluster), *(source_cluster)));
    }
};
#endif
