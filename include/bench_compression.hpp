#ifndef HTOOL_BENCHMARKS_BENCH_COMPRESSION_HPP
#define HTOOL_BENCHMARKS_BENCH_COMPRESSION_HPP

#include "geometry.hpp"
#include "kernel.hpp"
#include "misc.hpp"
#include <tuple>

template <typename T, template <typename> class Compressor>
std::tuple<double, double, double, double, int> bench_compression(int nr, int nc, double k, const std::vector<double> &p1, const std::vector<double> &p2, htool::VirtualCluster *cluster_target, htool::VirtualCluster *cluster_source, int vectorisation) {
    std::cout << nr << " " << nc << std::endl;

    // Htool parameters
    double epsilon = 1e-5;

    // Clustering
    std::vector<int> permutation_target, permutation_source;
    permutation_target = cluster_target->get_perm();
    permutation_source = cluster_source->get_perm();
    std::vector<double> p1_copy(p1), p2_copy(p2);
    if (vectorisation == 1 || vectorisation == 2) {
        std::vector<double> tmp(3 * nr, 0);
        for (int i = 0; i < nr; i++) {
            tmp[i + 0 * nr] = p1[permutation_target[i] * 3 + 0];
            tmp[i + 1 * nr] = p1[permutation_target[i] * 3 + 1];
            tmp[i + 2 * nr] = p1[permutation_target[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), tmp.size(), p1_copy.begin());
        tmp.resize(3 * nc);
        for (int i = 0; i < nc; i++) {
            tmp[i + 0 * nc] = p2[permutation_source[i] * 3 + 0];
            tmp[i + 1 * nc] = p2[permutation_source[i] * 3 + 1];
            tmp[i + 2 * nc] = p2[permutation_source[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), tmp.size(), p2_copy.begin());
    } else if (vectorisation == 3) {
        constexpr std::size_t simd_size = xsimd::simd_type<double>::size;
        int vec_size_target             = nr + (nr % simd_size);
        int vec_size_source             = nc + (nc % simd_size);
        p1_copy.resize(3 * vec_size_target, 0);
        p2_copy.resize(3 * vec_size_source, 0);
        std::vector<double> tmp(3 * nr, 0);
        for (int i = 0; i < nr; i++) {
            tmp[i + 0 * nr] = p1[permutation_target[i] * 3 + 0];
            tmp[i + 1 * nr] = p1[permutation_target[i] * 3 + 1];
            tmp[i + 2 * nr] = p1[permutation_target[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), nr, p1_copy.begin());
        std::copy_n(tmp.begin() + nr, nr, p1_copy.begin() + vec_size_target);
        std::copy_n(tmp.begin() + 2 * nr, nr, p1_copy.begin() + 2 * vec_size_target);
        tmp.resize(3 * nc);
        for (int i = 0; i < nc; i++) {
            tmp[i + 0 * nc] = p2[permutation_source[i] * 3 + 0];
            tmp[i + 1 * nc] = p2[permutation_source[i] * 3 + 1];
            tmp[i + 2 * nc] = p2[permutation_source[i] * 3 + 2];
        }
        std::copy_n(tmp.begin(), nc, p2_copy.begin());
        std::copy_n(tmp.begin() + nc, nc, p2_copy.begin() + vec_size_source);
        std::copy_n(tmp.begin() + 2 * nc, nc, p2_copy.begin() + 2 * vec_size_source);
    }

    // Generator
    std::unique_ptr<MyVirtualGenerator<T>> A;
    if (vectorisation == 0) {
        std::vector<double> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 0>>(3, nr, nc, p1_allocated, p2_allocated, k);
    } else if (vectorisation == 1) {
        std::vector<double> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 1>>(3, nr, nc, p1_allocated, p2_allocated, k);
    } else if (vectorisation == 2) {
        std::vector<double> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 2>>(3, nr, nc, p1_allocated, p2_allocated, k);
    } else if (vectorisation == 3) {
        constexpr std::size_t simd_size = xsimd::simd_type<double>::size;
        int vec_size_target             = nr + (nr % simd_size);
        int vec_size_source             = nc + (nc % simd_size);
        std::vector<double, xsimd::aligned_allocator<double>> p1_allocated(p1_copy.begin(), p1_copy.end()), p2_allocated(p2_copy.begin(), p2_copy.end());
        A = std::make_unique<MyGenerator<T, 3>>(3, nr, nc, p1_allocated, p2_allocated, k, vec_size_target, vec_size_source);
    }

    htool::underlying_type<T> norm_A = A->normFrob();

    // Timing
    double build_time, mat_prod_time, time(MPI_Wtime());

    // Compression
    time = MPI_Wtime();
    std::vector<int> target_coord(permutation_target), source_coord(permutation_source);
    if (vectorisation > 0) {
        std::iota(target_coord.begin(), target_coord.end(), int(0));
        std::iota(source_coord.begin(), source_coord.end(), int(0));
    }
    htool::LowRankMatrix<T> A_compressed(A->get_dimension(), target_coord, source_coord, -1, epsilon);

    std::shared_ptr<Compressor<T>> compressor = std::make_shared<Compressor<T>>();
    A_compressed.build(*A, *compressor, *cluster_target, p1_copy.data(), *cluster_target, p2_copy.data());

    build_time = MPI_Wtime() - time;
    std::cout << "Time  partialACA: " << build_time << std::endl;
    std::cout << "Error partialACA: " << std::setprecision(10) << Frobenius_absolute_error(A_compressed, *A) / norm_A << std::endl;
    std::cout << "Rank  partialACA: " << A_compressed.rank_of() << std::endl;

    std::vector<T> ones(nc, 1);
    time                         = MPI_Wtime();
    std::vector<T> out           = A_compressed * ones;
    mat_prod_time                = MPI_Wtime() - time;
    std::vector<int> perm_target = cluster_target->get_perm();
    std::vector<T> out_perm(out.size());
    for (int i = 0; i < perm_target.size(); i++) {
        out_perm[perm_target[i]] = out[i];
    }

    int line_number = (vectorisation > 0) ? std::distance(permutation_target.begin(), std::find(permutation_target.begin(), permutation_target.end(), 41)) : 41;
    std::cout << std::setprecision(10);
    std::cout << out_perm[41] / T(nc) << " " << A->sanity_check(line_number) << std::endl;

    std::tuple<double, double, double, double, int> results(build_time, mat_prod_time, 0., 0., A_compressed.rank_of());
    if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
        std::get<2>(results) = out_perm[41] / T(nc);
    } else if constexpr (std::is_same_v<T, std::complex<float>> || std::is_same_v<T, std::complex<double>>) {
        std::get<2>(results) = std::real(out_perm[41] / T(nc));
        std::get<3>(results) = std::imag(out_perm[41] / T(nc));
    }

    return results;
}

#endif