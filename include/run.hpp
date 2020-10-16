#ifndef HTOOL_BENCHMARKS_RUN
#define HTOOL_BENCHMARKS_RUN

#include "geometry.hpp"
#include "kernels.hpp"
#include <htool/htool.hpp>

template <typename T, template <typename, typename> class LowRankMatrix, class ClusterImpl, template <typename> class AdmissibleCondition, typename Kernel>
void run(int argc, char *argv[]) {

    if (argc != 2 && argc != 3) {
        std::cout << "usage: " << argv[0] << " <filename of geometry>\n";
        std::cout << "or\n";
        std::cout << "usage: " << argv[0] << " <filename of target geometry> <filename of source geometry>\n";
    } else {

        // Initialize the MPI environment
        MPI_Init(&argc, &argv);

        // Load geometry
        std::string filename_target;
        std::string filename_source;
        filename_target = std::string(argv[1]);
        filename_source = filename_target;
        if (argc == 3) {
            filename_source = std::string(argv[1]);
        }
        std::vector<std::array<double, 3>> geometry_target = LoadGeometry(filename_target);
        std::vector<std::array<double, 3>> geometry_source = LoadGeometry(filename_source);

        // Build kernel
        Kernel kernel(geometry_target, geometry_source);

        // Build HMatrix
        htool::HMatrix<T, LowRankMatrix, ClusterImpl, AdmissibleCondition> hmatrix(kernel, geometry_target, geometry_source);

        hmatrix.print_infos();

        // Finalize the MPI environment.
        MPI_Finalize();
    }
}

#endif