#ifndef CLI_HPP
#define CLI_HPP

#include "htool/clustering/implementations/partitioning.hpp"
#include "htool/clustering/interfaces/virtual_partitioning.hpp"
#include "htool/hmatrix/interfaces/virtual_generator.hpp"
#include "htool/hmatrix/interfaces/virtual_lrmat_generator.hpp"
#include "htool/hmatrix/lrmat/SVD.hpp"
#include "htool/hmatrix/lrmat/fullACA.hpp"
#include "htool/hmatrix/lrmat/partialACA.hpp"
#include "htool/hmatrix/lrmat/sympartialACA.hpp"
#include <cstddef>
#include <iostream>
#include <memory>

inline bool check_input_size(int argc, char *argv[]) {
    // Check for valid number of arguments
    if (argc < 1 || argc > 7) {
        std::cerr << "Usage: " << argv[0] << " <test_case> = {pbl_size|thread|ratio}[default=pbl_size] <symmetry_type>[default=N] = {N|S|H} <generator_type>[default=Laplace] = {Laplace|Helmholtz} <clustering_type>[default=PCARegularN] = {BoundingBoxRegular1|BoundingBoxGeometric1|PCARegular1|PCABoxGeometric1|BoundingBoxRegularN|BoundingBoxGeometricN|PCARegularN|PCABoxGeometricN} <low_rank_generator_type>[default=sympartialACA] = {SVD|fullACA|partialACA|sympartialACA} <policy>[default=par] = {seq|par|omp_task}" << std::endl;
        return true;
    }
    return false;
}

inline bool check_test_case_type(int argc, char *argv[], std::string &test_case_type) {
    test_case_type = (argc == 2) ? argv[1] : "pbl_size";
    if (test_case_type != "pbl_size" && test_case_type != "thread" && test_case_type != "ratio") {
        std::cerr << "Error: invalid test case type. Must be either \"pbl_size\", \"thread\" or \"ratio\"." << std::endl;
        return true;
    }
    return false;
}

inline bool check_symmetry_type(int argc, char *argv[], char &symmetry_type) {
    symmetry_type = (argc == 3) ? std::toupper(argv[2][0]) : 'N';
    if (symmetry_type != 'N' && symmetry_type != 'S') {
        std::cerr << "Error: invalid symmetry type. Must be either 'N' or 'S'." << std::endl;
        return true;
    }
    return false;
}

inline bool check_generator_type(int argc, char *argv[], std::string &generator_type) {
    generator_type = (argc == 4) ? argv[3] : "Laplace";
    if (generator_type != "Laplace" && generator_type != "Helmholtz") {
        std::cerr << "Error: invalid generator type. Must be either \"Laplace\" or \"Helmholtz\"." << std::endl;
        return true;
    }
    return false;
}

inline bool check_clustering_type(int argc, char *argv[], std::string &clustering_type) {
    clustering_type = (argc == 5) ? argv[4] : "PCARegularN";
    if (clustering_type != "BoundingBoxRegular1" && clustering_type != "BoundingBoxGeometric1" && clustering_type != "PCARegular1" && clustering_type != "PCABoxGeometric1" && clustering_type != "BoundingBoxRegularN" && clustering_type != "BoundingBoxGeometricN" && clustering_type != "PCARegularN" && clustering_type != "PCABoxGeometricN") {
        std::cerr << "Error: invalid clustering type. Must be either \"BoundingBoxRegular1\", \"BoundingBoxGeometric1\", \"PCARegular1\", \"PCABoxGeometric1\", \"BoundingBoxRegularN\", \"BoundingBoxGeometricN\", \"PCARegularN\" or \"PCABoxGeometricN\"." << std::endl;
        return true;
    }
    return false;
}

inline bool check_low_rank_generator_type(int argc, char *argv[], std::string &low_rank_generator_type) {
    low_rank_generator_type = (argc == 6) ? argv[5] : "sympartialACA";
    if (low_rank_generator_type != "sympartialACA" && low_rank_generator_type != "SVD" && low_rank_generator_type != "fullACA" && low_rank_generator_type != "partialACA" && low_rank_generator_type != "sympartialACA") {
        std::cerr << "Error: invalid clustering type. Must be either \"sympartialACA\", \"SVD\", \"fullACA\", \"partialACA\", or \"sympartialACA\"." << std::endl;
        return true;
    }
    return false;
}

inline bool check_policy_type(int argc, char *argv[], std::string &policy_type) {
    policy_type = (argc == 7) ? argv[6] : "par";
    if (policy_type != "seq" && policy_type != "par" && policy_type != "omp_task") {
        std::cerr << "Error: invalid policy type. Must be either \"seq\", \"par\", or \"task\"." << std::endl;
        return true;
    }
    return false;
}

template <typename CoordinatePrecision>
std::shared_ptr<htool::VirtualPartitioning<CoordinatePrecision>> process_clustering_type(std::string clustering_type) {
    if (clustering_type == "BoundingBoxRegular1") {
        return std::make_shared<htool::Partitioning<CoordinatePrecision, htool::ComputeBoundingBox<CoordinatePrecision>, htool::RegularSplitting<CoordinatePrecision>>>();
    } else if (clustering_type == "BoundingBoxGeometric1") {
        return std::make_shared<htool::Partitioning<CoordinatePrecision, htool::ComputeBoundingBox<CoordinatePrecision>, htool::GeometricSplitting<CoordinatePrecision>>>();
    } else if (clustering_type == "PCARegular1") {
        return std::make_shared<htool::Partitioning<CoordinatePrecision, htool::ComputeLargestExtent<CoordinatePrecision>, htool::RegularSplitting<CoordinatePrecision>>>();
    } else if (clustering_type == "PCABoxGeometric1") {
        return std::make_shared<htool::Partitioning<CoordinatePrecision, htool::ComputeLargestExtent<CoordinatePrecision>, htool::GeometricSplitting<CoordinatePrecision>>>();
    } else if (clustering_type == "BoundingBoxRegularN") {
        return std::make_shared<htool::Partitioning_N<CoordinatePrecision, htool::ComputeBoundingBox<CoordinatePrecision>, htool::RegularSplitting<CoordinatePrecision>>>();
    } else if (clustering_type == "BoundingBoxGeometricN") {
        return std::make_shared<htool::Partitioning_N<CoordinatePrecision, htool::ComputeBoundingBox<CoordinatePrecision>, htool::GeometricSplitting<CoordinatePrecision>>>();
    } else if (clustering_type == "PCARegularN") {
        return std::make_shared<htool::Partitioning_N<CoordinatePrecision, htool::ComputeLargestExtent<CoordinatePrecision>, htool::RegularSplitting<CoordinatePrecision>>>();
    } else if (clustering_type == "PCABoxGeometricN") {
        return std::make_shared<htool::Partitioning_N<CoordinatePrecision, htool::ComputeLargestExtent<CoordinatePrecision>, htool::GeometricSplitting<CoordinatePrecision>>>();
    } else {
    }
    return nullptr;
}

template <typename CoefficientPrecision>
std::shared_ptr<htool::VirtualInternalLowRankGenerator<CoefficientPrecision>> process_low_rank_generator_type(std::string low_rank_generator_type, htool::VirtualGenerator<CoefficientPrecision> &generator, const int *target_permutation, const int *source_permutation) {
    if (low_rank_generator_type == "SVD") {
        return std::make_shared<htool::SVD<CoefficientPrecision>>(generator, target_permutation, source_permutation);
    } else if (low_rank_generator_type == "fullACA") {
        return std::make_shared<htool::fullACA<CoefficientPrecision>>(generator, target_permutation, source_permutation);
    } else if (low_rank_generator_type == "partialACA") {
        return std::make_shared<htool::partialACA<CoefficientPrecision>>(generator, target_permutation, source_permutation);
    } else if (low_rank_generator_type == "sympartialACA") {
        return std::make_shared<htool::sympartialACA<CoefficientPrecision>>(generator, target_permutation, source_permutation);
    }
    return nullptr;
}

template <typename CoefficientPrecision>
std::shared_ptr<htool::VirtualInternalLowRankGenerator<CoefficientPrecision>> process_low_rank_generator_type(std::string low_rank_generator_type, htool::VirtualInternalGenerator<CoefficientPrecision> &generator) {
    if (low_rank_generator_type == "SVD") {
        return std::make_shared<htool::SVD<CoefficientPrecision>>(generator);
    } else if (low_rank_generator_type == "fullACA") {
        return std::make_shared<htool::fullACA<CoefficientPrecision>>(generator);
    } else if (low_rank_generator_type == "partialACA") {
        return std::make_shared<htool::partialACA<CoefficientPrecision>>(generator);
    } else if (low_rank_generator_type == "sympartialACA") {
        return std::make_shared<htool::sympartialACA<CoefficientPrecision>>(generator);
    }
    return nullptr;
}

#endif
