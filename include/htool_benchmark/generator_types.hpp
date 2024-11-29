#ifndef HTOOL_BENCHMARK_GENERATOR_TYPE_HPP
#define HTOOL_BENCHMARK_GENERATOR_TYPE_HPP
#include <cmath>
#include <htool/hmatrix/interfaces/virtual_generator.hpp>
#include <numeric>
#include <vector>
namespace htool_benchmark {

class LaplaceLikeGenerator : public htool::VirtualGenerator<double> {
    int m_space_dimension;
    const std::vector<double> &m_target_points;
    const std::vector<double> &m_source_points;

  public:
    LaplaceLikeGenerator(int dimension, const std::vector<double> &target_points, const std::vector<double> &source_points) : m_space_dimension(dimension), m_target_points(target_points), m_source_points(source_points) {}

    void copy_submatrix(int M, int N, const int *rows, const int *cols, double *ptr) const override {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                ptr[i + M * j] = this->get_coef(rows[i], cols[j]);
            }
        }
    }

    double get_coef(const int &i, const int &j) const {
        return 1. / (1e-5 + std::sqrt(std::inner_product(m_target_points.begin() + m_space_dimension * i, m_target_points.begin() + m_space_dimension * i + m_space_dimension, m_source_points.begin() + m_space_dimension * j, double(0), std::plus<double>(), [](double u, double v) { return (u - v) * (u - v); })));
    }
};

class HelmholtzLikeGenerator : public htool::VirtualGenerator<std::complex<double>> { // TODO : valider
    int m_space_dimension;
    const std::vector<double> &m_target_points;
    const std::vector<double> &m_source_points;

  public:
    HelmholtzLikeGenerator(int dimension, const std::vector<double> &target_points, const std::vector<double> &source_points) : m_space_dimension(dimension), m_target_points(target_points), m_source_points(source_points) {}

    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, std::complex<double> *ptr) const override {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                ptr[i + M * j] = this->get_coef(rows[i], cols[j]);
            }
        }
    }

    std::complex<double> get_coef(const int &i, const int &j) const {
        double n_xMy = std::sqrt(std::inner_product(m_target_points.begin() + m_space_dimension * i, m_target_points.begin() + m_space_dimension * i + m_space_dimension, m_source_points.begin() + m_space_dimension * j, double(0), std::plus<double>(), [](double u, double v) { return (u - v) * (u - v); }));
        return (cos(10 * n_xMy) + std::complex<double>(0, 1) * sin(10 * n_xMy)) / (1e-5 + n_xMy);
    }
};

} // namespace htool_benchmark

#endif
