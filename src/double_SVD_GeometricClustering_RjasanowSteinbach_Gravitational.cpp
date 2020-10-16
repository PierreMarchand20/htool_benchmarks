#include "kernels.hpp"
#include "run.hpp"

int main(int argc, char *argv[]) {

    run<double, htool::SVD, htool::GeometricClustering, htool::RjasanowSteinbach, GravitationalKernel>(argc, argv);

    return 0;
}
