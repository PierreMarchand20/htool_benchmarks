#ifndef HTOOL_BENCHMARKS_KERNELS
#define HTOOL_BENCHMARKS_KERNELS

#include <htool/htool.hpp>

typedef std::array<double, 3> R3;

class GravitationalKernel : public htool::IMatrix<double> {
    const std::vector<R3> &p1;
    const std::vector<R3> &p2;

  public:
    GravitationalKernel(const std::vector<R3> &p10, const std::vector<R3> &p20) : IMatrix<double>(p10.size(), p20.size()), p1(p10), p2(p20) {}

    double get_coef(const int &i, const int &j) const { return 1. / ((p1[i][0] - p1[j][0]) * (p1[i][0] - p1[j][0]) + (p1[i][1] - p1[j][1]) * (p1[i][1] - p1[j][1]) + (p1[i][2] - p1[j][2]) * (p1[i][2] - p1[j][2])); }

    std::vector<double> operator*(std::vector<double> a) {
        std::vector<double> result(p1.size(), 0);
        for (int i = 0; i < p1.size(); i++) {
            for (int k = 0; k < p2.size(); k++) {
                result[i] += this->get_coef(i, k) * a[k];
            }
        }
        return result;
    }

    double normFrob() {
        double norm = 0;
        for (int i = 0; i < p1.size(); i++) {
            for (int k = 0; k < p2.size(); k++) {
                norm = norm + std::pow(this->get_coef(i, k), 2);
            }
        }
        return sqrt(norm);
    }
};

#endif