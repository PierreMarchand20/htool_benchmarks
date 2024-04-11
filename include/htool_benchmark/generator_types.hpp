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

class OptimizedLaplaceLikeGenerator : public htool::VirtualInternalGenerator<double> {
    int m_space_dimension;
    const std::vector<double> &m_target_points;
    const std::vector<double> &m_source_points;

  public:
    OptimizedLaplaceLikeGenerator(int dimension, const std::vector<double> &target_points, const std::vector<double> &source_points) : m_space_dimension(dimension), m_target_points(target_points), m_source_points(source_points) {}

    void copy_submatrix(int M, int N, int row_offset, int col_offset, double *ptr) const override {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < M; i++) {
                std::array<double, 3> distance;
                for (int k = 0; k < 3; k++) {
                    distance[k] = m_target_points[(i + row_offset) + m_target_points.size() / 3 * k] - m_source_points[(j + col_offset) + m_source_points.size() / 3 * k];
                }
                double squared_distance{0};
                for (int k = 0; k < 3; k++) {
                    squared_distance += distance[k] * distance[k];
                }
                ptr[i + M * j] = 1. / (sqrt(squared_distance) + 1e-5);
            }
        }
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

class OptimizedHelmholtzLikeGenerator : public htool::VirtualInternalGenerator<std::complex<double>> {
    int m_space_dimension;
    const std::vector<double> &m_target_points;
    const std::vector<double> &m_source_points;

  public:
    OptimizedHelmholtzLikeGenerator(int dimension, const std::vector<double> &target_points, const std::vector<double> &source_points) : m_space_dimension(dimension), m_target_points(target_points), m_source_points(source_points) {}

    void copy_submatrix(int M, int N, int row_offset, int col_offset, std::complex<double> *ptr) const override {
        int nb_target_points = m_target_points.size() / 3;
        int nb_source_points = m_source_points.size() / 3;
        double dx, dy, dz;
        double squared_distance, distance, real_part, imag_part;
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < M; i++) {
                dx = m_target_points[(i + row_offset) + nb_target_points * 0] - m_source_points[(j + col_offset) + nb_source_points * 0];
                dy = m_target_points[(i + row_offset) + nb_target_points * 1] - m_source_points[(j + col_offset) + nb_source_points * 1];
                dz = m_target_points[(i + row_offset) + nb_target_points * 2] - m_source_points[(j + col_offset) + nb_source_points * 2];
                squared_distance = dx*dx+dy*dy+dz*dz;
                distance = sqrt(squared_distance);
                real_part = cos(10. * distance) / (distance+1e-5);
                imag_part = sin(10. * distance) / (distance+1e-5);
                ptr[i + M * j] = std::complex<double>(real_part,imag_part);
            }
        }
        // for (int j = 0; j < N; j++) {
        //     for (int i = 0; i < M; i++) {
        //         std::array<double, 3> distances;
        //         for (int k = 0; k < 3; k++) {
        //             distances[k] = m_target_points[(i + row_offset) + m_target_points.size() / 3 * k] - m_source_points[(j + col_offset) + m_source_points.size() / 3 * k];
        //         }
        //         double squared_distance{0};
        //         for (int k = 0; k < 3; k++) {
        //             squared_distance += distances[k] * distances[k];
        //         }
        //         // double distance  = sqrt(squared_distance);
        //         // double real_part = cos(10 * distance)/(distance+1e-5);
        //         // double imag_part = sin(10 * distance)/(distance+1e-5); 
        //         // ptr[i + M * j]  = real_part + imag_part;
        //     }
        // }
    }
};

} // namespace htool_benchmark

#endif
