#ifndef HTOOL_BENCHMARKS_BENCH_HMATRIX_HPP
#define HTOOL_BENCHMARKS_BENCH_HMATRIX_HPP

#include "geometry.hpp"
#include "kernel.hpp"
#include "misc.hpp"
#include <tuple>

template <typename T, template <typename> class Compressor>
std::tuple<double, double, double, double, double> bench_hmatrix(int n, double k, const std::vector<double> &p1, const std::vector<double> &p2, std::shared_ptr<htool::VirtualCluster> cluster_target, std::shared_ptr<htool::VirtualCluster> cluster_source, int vectorisation) {

    // Htool parameters
    double epsilon        = 1e-5;
    double eta            = 3.0;
    double minclustersize = 200;

    // Clustering
    std::vector<int> permutation_target, permutation_source;
    permutation_target = cluster_target->get_perm();
    std::vector<double> p1_copy(p1), p2_copy(p2);
    if (vectorisation == 1 || vectorisation == 2) {
        permutation_source = cluster_source->get_perm();
        std::vector<double> tmp(3 * n, 0);
        for (int i = 0; i < n; i++) {
            tmp[i + 0 * n] = p1[permutation_target[i] * 3 + 0];
            tmp[i + 1 * n] = p1[permutation_target[i] * 3 + 1];
            tmp[i + 2 * n] = p1[permutation_target[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), tmp.size(), p1_copy.begin());
        for (int i = 0; i < n; i++) {
            tmp[i + 0 * n] = p2[permutation_source[i] * 3 + 0];
            tmp[i + 1 * n] = p2[permutation_source[i] * 3 + 1];
            tmp[i + 2 * n] = p2[permutation_source[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), tmp.size(), p2_copy.begin());
    } else if (vectorisation == 3) {
        constexpr std::size_t simd_size = xsimd::simd_type<double>::size;

        int vec_size = n + (n % simd_size);
        p1_copy.resize(3 * vec_size, 0);
        p2_copy.resize(3 * vec_size, 0);
        permutation_source = cluster_source->get_perm();
        std::vector<double> tmp(3 * n, 0);
        for (int i = 0; i < n; i++) {
            tmp[i + 0 * n] = p1[permutation_target[i] * 3 + 0];
            tmp[i + 1 * n] = p1[permutation_target[i] * 3 + 1];
            tmp[i + 2 * n] = p1[permutation_target[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), n, p1_copy.begin());
        std::copy_n(tmp.begin() + n, n, p1_copy.begin() + vec_size);
        std::copy_n(tmp.begin() + 2 * n, n, p1_copy.begin() + 2 * vec_size);
        for (int i = 0; i < n; i++) {
            tmp[i + 0 * n] = p2[permutation_source[i] * 3 + 0];
            tmp[i + 1 * n] = p2[permutation_source[i] * 3 + 1];
            tmp[i + 2 * n] = p2[permutation_source[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), n, p2_copy.begin());
        std::copy_n(tmp.begin() + n, n, p2_copy.begin() + vec_size);
        std::copy_n(tmp.begin() + 2 * n, n, p2_copy.begin() + 2 * vec_size);
    }

    // Generator
    std::unique_ptr<MyVirtualGenerator<T>> A;
    if (vectorisation == 0) {
        std::vector<double> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 0>>(3, n, n, p1_allocated, p2_allocated, k);
    } else if (vectorisation == 1) {
        std::vector<double> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 1>>(3, n, n, p1_allocated, p2_allocated, k);
    } else if (vectorisation == 2) {
        std::vector<double> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 2>>(3, n, n, p1_allocated, p2_allocated, k);
    } else if (vectorisation == 3) {
        constexpr std::size_t simd_size = xsimd::simd_type<double>::size;
        int vec_size                    = n + (n % simd_size);
        std::vector<double, xsimd::aligned_allocator<double>> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 3>>(3, n, n, p1_allocated, p2_allocated, k, vec_size, vec_size);
    }

    htool::underlying_type<T> norm_A = A->normFrob();

    // Hmatrix
    std::vector<T> x(n, 1), result(n, 0);
    double time1 = MPI_Wtime();
    htool::HMatrix<T, Compressor, htool::RjasanowSteinbach> HA(cluster_target, cluster_source, epsilon, eta, 'N', 'N');
    if (vectorisation > 0) {
        HA.set_use_permutation(false);
    }
    HA.build_auto(*A, p1_copy.data(), p2_copy.data());
    time1 = MPI_Wtime() - time1;
    HA.print_infos();
    HA.set_use_permutation(false);
    double time2 = MPI_Wtime();
    for (int i = 0; i < 25; i++) {
        HA.mvprod_global_to_global(x.data(), result.data());
    }
    time2 = MPI_Wtime() - time2;

    std::cout << "Total building time: " << time1 << std::endl;
    std::cout << "Total matvec time: " << time2 << std::endl;

    std::vector<T> tmp(n, 0);
    for (int i = 0; i < n; i++) {
        tmp[permutation_target[i]] = result[i];
    }
    std::copy_n(tmp.begin(), tmp.size(), result.begin());

    int line_number = (vectorisation > 0) ? std::distance(permutation_target.begin(), std::find(permutation_target.begin(), permutation_target.end(), 41)) : 41;
    std::cout << std::setprecision(10);
    std::cout << result[41] / T(n) << " " << A->sanity_check(line_number) << std::endl;

    // Save
    std::tuple<double, double, double, double, double> results(time1, time2, 0., 0., HA.space_saving());

    if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
        std::get<2>(results) = result[41] / T(n);
    } else if constexpr (std::is_same_v<T, std::complex<float>> || std::is_same_v<T, std::complex<double>>) {
        std::get<2>(results) = std::real(result[41] / T(n));
        std::get<3>(results) = std::imag(result[41] / T(n));
    }

    return results;
}

#endif