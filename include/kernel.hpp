#ifndef HTOOL_BENCHMARKS_KERNEL_HPP
#define HTOOL_BENCHMARKS_KERNEL_HPP

#include <complex>
#include <htool/htool.hpp>
#include <iomanip>
#include <omp.h>

#ifndef __INTEL_COMPILER
#    include <vectorclass.h>
#    include "vectormath_exp.h"
#    include "vectormath_trig.h"
#    include <complexvec1.h>
#endif
#include "misc.hpp"
#include "xsimd/xsimd.hpp"


template <typename T>
class MyVirtualGenerator : public htool::VirtualGenerator<T> {
  protected:
    htool::underlying_type<T> freq;

  public:
    // Constructor
    MyVirtualGenerator(int space_dim0, int nr, int nc, htool::underlying_type<T> freq0)
        : htool::VirtualGenerator<T>(nr, nc), freq(freq0) {}

    // Virtual function to overload
    T get_coef(const int &k, const int &j) const {}

    // Virtual function to overload
    virtual void copy_submatrix(int M, int N, const int *const rows, const int *const cols, T *ptr) const = 0;

    T sanity_check(int line_number) {
        std::vector<T> line(this->nc);
        std::vector<int> col_numbers(this->nc);
        std::iota(col_numbers.begin(), col_numbers.end(), int(0));
        this->copy_submatrix(1, this->nc, &line_number, col_numbers.data(), line.data());
        return std::accumulate(line.begin(), line.end(), T(0)) / T(this->nc);
    }

    htool::underlying_type<T> normFrob() {
        std::vector<T> mat(this->nr * this->nc);
        std::vector<int> line_numbers(this->nr), col_numbers(this->nc);
        std::iota(col_numbers.begin(), col_numbers.end(), int(0));
        std::iota(line_numbers.begin(), line_numbers.end(), int(0));
        this->copy_submatrix(this->nr, this->nc, line_numbers.data(), col_numbers.data(), mat.data());
        return htool::norm2(mat);
    }

    htool::underlying_type<T> get_freq() const { return this->freq; };
};

template <typename T, int vectorized>
class MyGenerator : public MyVirtualGenerator<T> {
    std::vector<double, allocator_type<double, vectorized>> points_target;
    std::vector<double, allocator_type<double, vectorized>> points_source;
    int space_dim;
    int inct;
    int incs;

  public:
    // Constructor
    MyGenerator(int space_dim0, int nr, int nc, const std::vector<double, allocator_type<double, vectorized>> &p10, const std::vector<double, allocator_type<double, vectorized>> &p20, htool::underlying_type<T> freq0, std::size_t inct0 = 0, std::size_t incs0 = 0)
        : MyVirtualGenerator<T>(space_dim0, nr, nc, freq0), points_target(p10), points_source(p20), space_dim(space_dim0), inct(inct0), incs(incs0) {
    }

    // Virtual function to overload
    T get_coef(const int &k, const int &j) const {}

    // Virtual function to overload
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, T *ptr) const {};
};

// No vectorization

template <>
void MyGenerator<float, 0>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, float *ptr) const {
    double dx, dy, dz;
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++) {
            dx             = this->points_target[this->space_dim * rows[j] + 0] - this->points_source[this->space_dim * cols[k] + 0];
            dy             = this->points_target[this->space_dim * rows[j] + 1] - this->points_source[this->space_dim * cols[k] + 1];
            dz             = this->points_target[this->space_dim * rows[j] + 2] - this->points_source[this->space_dim * cols[k] + 2];
            ptr[j + k * M] = 1.f / ((4 * 3.14159265358979f) * (sqrt(dx * dx + dy * dy + dz * dz)));
        }
    }
}

template <>
void MyGenerator<double, 0>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const {
    double dx, dy, dz;
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++) {
            dx             = this->points_target[space_dim * rows[j] + 0] - this->points_source[space_dim * cols[k] + 0];
            dy             = this->points_target[space_dim * rows[j] + 1] - this->points_source[space_dim * cols[k] + 1];
            dz             = this->points_target[space_dim * rows[j] + 2] - this->points_source[space_dim * cols[k] + 2];
            ptr[j + k * M] = 1. / ((4 * 3.141592653589793238463) * (sqrt(dx * dx + dy * dy + dz * dz)));
        }
    }
}

template <>
void MyGenerator<std::complex<float>, 0>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, std::complex<float> *ptr) const {
    double dx, dy, dz;
    float d;
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++) {
            dx             = this->points_target[space_dim * rows[j] + 0] - this->points_source[space_dim * cols[k] + 0];
            dy             = this->points_target[space_dim * rows[j] + 1] - this->points_source[space_dim * cols[k] + 1];
            dz             = this->points_target[space_dim * rows[j] + 2] - this->points_source[space_dim * cols[k] + 2];
            d              = sqrt(dx * dx + dy * dy + dz * dz);
            ptr[j + k * M] = (1.f / (4 * 3.14159265358979f)) * exp(std::complex<float>(0, 1) * freq * d) / (d);
        }
    }
}

template <>
void MyGenerator<std::complex<double>, 0>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, std::complex<double> *ptr) const {
    double dx, dy, dz, d;
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++) {
            dx             = this->points_target[space_dim * rows[j] + 0] - this->points_source[space_dim * cols[k] + 0];
            dy             = this->points_target[space_dim * rows[j] + 1] - this->points_source[space_dim * cols[k] + 1];
            dz             = this->points_target[space_dim * rows[j] + 2] - this->points_source[space_dim * cols[k] + 2];
            d              = sqrt(dx * dx + dy * dy + dz * dz);
            ptr[j + k * M] = (1.f / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
        }
    }
}

// vector-class-library Agner Fog's library
#ifndef __INTEL_COMPILER
template <>
void MyGenerator<float, 1>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, float *ptr) const {
    double ddx, ddy, ddz;
    Vec8d charge(1. / (4 * 3.141592653589793238463));
    Vec8d Xt_vec, Yt_vec, Zt_vec;
    Vec8f pot_vec(0);
    Vec8d dx, dy, dz, r2;
    if (M >= 8) {
        std::size_t vec_size = M - M % 8;
        for (std::size_t j = 0; j < vec_size; j += 8) {

            // load
            Xt_vec.load(this->points_target.data() + rows[0] + 0 * this->nr + j);
            Yt_vec.load(this->points_target.data() + rows[0] + 1 * this->nr + j);
            Zt_vec.load(this->points_target.data() + rows[0] + 2 * this->nr + j);

            for (std::size_t k = 0; k < N; k++) {
                dx      = Vec8d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec8d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec8d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = to_float(charge / sqrt(r2));
                pot_vec.store(&(ptr[j + k * M]));
            }
        }
        if (vec_size < M) {
            int padding = M - vec_size;
            // load
            Xt_vec.load_partial(padding, this->points_target.data() + rows[0] + 0 * this->nr + vec_size);
            Yt_vec.load_partial(padding, this->points_target.data() + rows[0] + 1 * this->nr + vec_size);
            Zt_vec.load_partial(padding, this->points_target.data() + rows[0] + 2 * this->nr + vec_size);
            for (std::size_t k = 0; k < N; k++) {
                dx      = Vec8d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec8d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec8d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = to_float(charge / sqrt(r2));
                pot_vec.store_partial(padding, &(ptr[vec_size + k * M]));
            }
        }
    } else if (N >= 8) {
        std::size_t vec_size = N - N % 8;
        std::vector<float> tmp(8);
        for (std::size_t k = 0; k < vec_size; k += 8) {
            // load
            Xt_vec.load(this->points_source.data() + cols[0] + 0 * this->nc + k);
            Yt_vec.load(this->points_source.data() + cols[0] + 1 * this->nc + k);
            Zt_vec.load(this->points_source.data() + cols[0] + 2 * this->nc + k);

            for (std::size_t j = 0; j < M; j++) {
                dx      = Vec8d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec8d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec8d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = to_float(charge / sqrt(r2));
                pot_vec.store(tmp.data());
                for (int i = 0; i < 8; i++) {
                    ptr[j + (k + i) * M] = tmp[i];
                }
            }
        }
        if (vec_size < N) {
            int padding = N - vec_size;
            // load
            Xt_vec.load_partial(padding, this->points_source.data() + cols[0] + 0 * this->nc + vec_size);
            Yt_vec.load_partial(padding, this->points_source.data() + cols[0] + 1 * this->nc + vec_size);
            Zt_vec.load_partial(padding, this->points_source.data() + cols[0] + 2 * this->nc + vec_size);
            for (std::size_t j = 0; j < M; j++) {
                dx      = Vec8d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec8d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec8d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = to_float(charge / sqrt(r2));
                pot_vec.store(tmp.data());
                for (std::size_t i = 0; i < padding; i++) {
                    ptr[j + (vec_size + i) * M] = tmp[i];
                }
            }
        }
    } else {
        for (std::size_t k = 0; k < N; k++) {
            for (std::size_t j = 0; j < M; j++) {
                ddx            = this->points_target[this->space_dim * rows[j] + 0] - this->points_source[this->space_dim * cols[k] + 0];
                ddy            = this->points_target[this->space_dim * rows[j] + 1] - this->points_source[this->space_dim * cols[k] + 1];
                ddz            = this->points_target[this->space_dim * rows[j] + 2] - this->points_source[this->space_dim * cols[k] + 2];
                ptr[j + k * M] = 1.f / ((4 * 3.14159265358979f) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    }
}

template <>
void MyGenerator<double, 1>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const {
    double ddx, ddy, ddz;
    Vec8d charge(1. / (4 * 3.141592653589793238463));
    Vec8d Xt_vec, Yt_vec, Zt_vec, pot_vec, r2;
    Vec8d dx, dy, dz;

    if (M >= 8) {
        std::size_t vec_size = M - M % 8;
        for (std::size_t j = 0; j < vec_size; j += 8) {

            // load
            Xt_vec.load(this->points_target.data() + rows[0] + 0 * this->nr + j);
            Yt_vec.load(this->points_target.data() + rows[0] + 1 * this->nr + j);
            Zt_vec.load(this->points_target.data() + rows[0] + 2 * this->nr + j);

            for (std::size_t k = 0; k < N; k++) {
                dx      = Vec8d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec8d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec8d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                pot_vec.store(&(ptr[j + k * M]));
            }
        }
        if (vec_size < M) {
            int padding = M - vec_size;

            // load
            Xt_vec.load_partial(padding, this->points_target.data() + rows[0] + 0 * this->nr + vec_size);
            Yt_vec.load_partial(padding, this->points_target.data() + rows[0] + 1 * this->nr + vec_size);
            Zt_vec.load_partial(padding, this->points_target.data() + rows[0] + 2 * this->nr + vec_size);
            for (std::size_t k = 0; k < N; k++) {
                dx      = Vec8d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec8d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec8d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                pot_vec.store_partial(padding, &(ptr[vec_size + k * M]));
            }
        }
    } else if (N >= 8) {
        std::vector<double> tmp(8);
        std::size_t vec_size = N - N % 8;
        for (std::size_t k = 0; k < vec_size; k += 8) {

            // load
            Xt_vec.load(this->points_source.data() + cols[0] + 0 * this->nc + k);
            Yt_vec.load(this->points_source.data() + cols[0] + 1 * this->nc + k);
            Zt_vec.load(this->points_source.data() + cols[0] + 2 * this->nc + k);

            for (std::size_t j = 0; j < M; j++) {
                dx      = Vec8d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec8d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec8d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                pot_vec.store(tmp.data());
                for (int i = 0; i < 8; i++) {
                    ptr[j + (k + i) * M] = tmp[i];
                }
            }
        }
        if (vec_size < N) {
            int padding = N - vec_size;

            // load
            Xt_vec.load_partial(padding, this->points_source.data() + cols[0] + 0 * this->nc + vec_size);
            Yt_vec.load_partial(padding, this->points_source.data() + cols[0] + 1 * this->nc + vec_size);
            Zt_vec.load_partial(padding, this->points_source.data() + cols[0] + 2 * this->nc + vec_size);
            for (std::size_t j = 0; j < M; j++) {
                dx      = Vec8d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec8d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec8d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                pot_vec.store(tmp.data());
                for (std::size_t i = 0; i < padding; i++) {
                    ptr[j + (vec_size + i) * M] = tmp[i];
                }
            }
        }
    } else {

        for (std::size_t k = 0; k < N; k++) {
            for (std::size_t j = 0; j < M; j++) {
                ddx            = this->points_target[this->space_dim * rows[j] + 0] - this->points_source[this->space_dim * cols[k] + 0];
                ddy            = this->points_target[this->space_dim * rows[j] + 1] - this->points_source[this->space_dim * cols[k] + 1];
                ddz            = this->points_target[this->space_dim * rows[j] + 2] - this->points_source[this->space_dim * cols[k] + 2];
                ptr[j + k * M] = 1 / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    }
}

template <>
void MyGenerator<std::complex<float>, 1>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, std::complex<float> *ptr) const {
    double ddx, ddy, ddz;
    Vec4d charge = Vec4d(1. / (4 * 3.141592653589793238463));
    Vec4d Xt_vec, Yt_vec, Zt_vec;
    Complex4f pot_vec;
    Vec4d dx, dy, dz, r, aux, cos, sin;
    float d;
    if (M >= 4) {
        std::size_t vec_size = M - M % 4;
        for (std::size_t j = 0; j < vec_size; j += 4) {

            // load
            Xt_vec.load(this->points_target.data() + rows[0] + 0 * this->nr + j);
            Yt_vec.load(this->points_target.data() + rows[0] + 1 * this->nr + j);
            Zt_vec.load(this->points_target.data() + rows[0] + 2 * this->nr + j);

            for (std::size_t k = 0; k < N; k++) {
                dx      = Vec4d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec4d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec4d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = to_float(Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin))));
                pot_vec.store(reinterpret_cast<float *>(&(ptr[j + k * M])));
            }
        }
        if (vec_size < M) {
            int padding = M - vec_size;

            // load
            Xt_vec.load_partial(padding, this->points_target.data() + rows[0] + 0 * this->nr + vec_size);
            Yt_vec.load_partial(padding, this->points_target.data() + rows[0] + 1 * this->nr + vec_size);
            Zt_vec.load_partial(padding, this->points_target.data() + rows[0] + 2 * this->nr + vec_size);
            for (int k = 0; k < N; k++) {
                dx      = Vec4d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec4d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec4d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = to_float(Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin))));
                pot_vec.to_vector().store_partial(2 * padding, reinterpret_cast<float *>(&(ptr[vec_size + k * M])));
            }
        }
    } else if (N >= 4) {
        std::size_t vec_size = N - N % 4;
        std::vector<std::complex<float>> tmp(4);
        for (std::size_t k = 0; k < vec_size; k += 4) {

            // load
            Xt_vec.load(this->points_source.data() + cols[0] + 0 * this->nc + k);
            Yt_vec.load(this->points_source.data() + cols[0] + 1 * this->nc + k);
            Zt_vec.load(this->points_source.data() + cols[0] + 2 * this->nc + k);

            for (std::size_t j = 0; j < M; j++) {
                dx      = Vec4d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec4d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec4d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = to_float(Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin))));
                pot_vec.store(reinterpret_cast<float *>(tmp.data()));
                for (int i = 0; i < 4; i++) {
                    ptr[j + (k + i) * M] = tmp[i];
                }
            }
        }
        if (vec_size < N) {
            int padding = N - vec_size;

            // load
            Xt_vec.load_partial(padding, this->points_source.data() + cols[0] + 0 * this->nc + vec_size);
            Yt_vec.load_partial(padding, this->points_source.data() + cols[0] + 1 * this->nc + vec_size);
            Zt_vec.load_partial(padding, this->points_source.data() + cols[0] + 2 * this->nc + vec_size);
            for (int j = 0; j < M; j++) {
                dx      = Vec4d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec4d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec4d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = to_float(Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin))));
                pot_vec.store(reinterpret_cast<float *>(tmp.data()));
                for (int i = 0; i < padding; i++) {
                    ptr[j + (vec_size + i) * M] = tmp[i];
                }
            }
        }
    } else {
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < M; j++) {
                dx             = this->points_target[space_dim * rows[j] + 0] - this->points_source[space_dim * cols[k] + 0];
                dy             = this->points_target[space_dim * rows[j] + 1] - this->points_source[space_dim * cols[k] + 1];
                dz             = this->points_target[space_dim * rows[j] + 2] - this->points_source[space_dim * cols[k] + 2];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1.f / (4 * 3.14159265358979f)) * exp(std::complex<float>(0, 1) * freq * d) / (d);
            }
        }
    }
}

template <>
void MyGenerator<std::complex<double>, 1>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, std::complex<double> *ptr) const {
    double ddx, ddy, ddz, d;
    Vec4d charge = Vec4d(1. / (4 * 3.141592653589793238463)), r, aux;
    Vec4d cos, sin;
    Vec4d Xt_vec, Yt_vec, Zt_vec;
    Complex4d pot_vec;
    Vec4d dx, dy, dz;

    if (M >= 4) {
        std::size_t vec_size = M - M % 4;
        for (std::size_t j = 0; j < vec_size; j += 4) {

            // load
            Xt_vec.load(this->points_target.data() + rows[0] + 0 * this->nr + j);
            Yt_vec.load(this->points_target.data() + rows[0] + 1 * this->nr + j);
            Zt_vec.load(this->points_target.data() + rows[0] + 2 * this->nr + j);

            for (std::size_t k = 0; k < N; k++) {
                dx      = Vec4d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec4d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec4d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin)));
                pot_vec.store(reinterpret_cast<double *>(&(ptr[j + k * M])));
            }
        }
        if (vec_size < M) {
            int padding = M - vec_size;

            // load
            Xt_vec.load_partial(padding, this->points_target.data() + rows[0] + 0 * this->nr + vec_size);
            Yt_vec.load_partial(padding, this->points_target.data() + rows[0] + 1 * this->nr + vec_size);
            Zt_vec.load_partial(padding, this->points_target.data() + rows[0] + 2 * this->nr + vec_size);
            for (int k = 0; k < N; k++) {
                dx      = Vec4d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = Vec4d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = Vec4d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin)));
                pot_vec.to_vector().store_partial(2 * padding, reinterpret_cast<double *>(&(ptr[vec_size + k * M])));
            }
        }
    } else if (N >= 4) {
        std::size_t vec_size = N - N % 4;
        std::vector<std::complex<double>> tmp(4);
        for (std::size_t k = 0; k < vec_size; k += 4) {

            // load
            Xt_vec.load(this->points_source.data() + cols[0] + 0 * this->nc + k);
            Yt_vec.load(this->points_source.data() + cols[0] + 1 * this->nc + k);
            Zt_vec.load(this->points_source.data() + cols[0] + 2 * this->nc + k);

            for (std::size_t j = 0; j < M; j++) {
                dx      = Vec4d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec4d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec4d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin)));
                pot_vec.store(reinterpret_cast<double *>(tmp.data()));
                for (int i = 0; i < 4; i++) {
                    ptr[j + (k + i) * M] = tmp[i];
                }
            }
        }
        if (vec_size < N) {
            int padding = N - vec_size;

            // load
            Xt_vec.load_partial(padding, this->points_source.data() + cols[0] + 0 * this->nc + vec_size);
            Yt_vec.load_partial(padding, this->points_source.data() + cols[0] + 1 * this->nc + vec_size);
            Zt_vec.load_partial(padding, this->points_source.data() + cols[0] + 2 * this->nc + vec_size);
            for (int j = 0; j < M; j++) {
                dx      = Vec4d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = Vec4d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = Vec4d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r       = sqrt(dx * dx + dy * dy + dz * dz);
                aux     = charge / r;
                sin     = sincos(&cos, freq * r);
                pot_vec = Complex4d(permute8<0, 4, 1, 5, 2, 6, 3, 7>(Vec8d(aux * cos, aux * sin)));
                pot_vec.store(reinterpret_cast<double *>(tmp.data()));
                for (int i = 0; i < padding; i++) {
                    ptr[j + (vec_size + i) * M] = tmp[i];
                }
            }
        }
    } else {
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < M; j++) {
                ddx            = this->points_target[space_dim * rows[j] + 0] - this->points_source[space_dim * cols[k] + 0];
                ddy            = this->points_target[space_dim * rows[j] + 1] - this->points_source[space_dim * cols[k] + 1];
                ddz            = this->points_target[space_dim * rows[j] + 2] - this->points_source[space_dim * cols[k] + 2];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1.f / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }
    }
}
#endif

// xsimd library (unaligned)

template <>
void MyGenerator<double, 2>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const {
    using b_type_d                  = xsimd::simd_type<double>;
    constexpr std::size_t simd_size = b_type_d::size;

    b_type_d charge(1. / (4 * 3.141592653589793238463)), r2;
    b_type_d dx, dy, dz;
    b_type_d Xt_vec, Yt_vec, Zt_vec;
    b_type_d pot_vec;
    double ddx, ddy, ddz;
    if (M >= simd_size) {

        std::size_t vec_size = M - M % simd_size;
        for (std::size_t j = 0; j < vec_size; j += simd_size) {

            // load
            Xt_vec.load_unaligned(this->points_target.data() + rows[0] + 0 * this->nr + j);
            Yt_vec.load_unaligned(this->points_target.data() + rows[0] + 1 * this->nr + j);
            Zt_vec.load_unaligned(this->points_target.data() + rows[0] + 2 * this->nr + j);

            for (std::size_t k = 0; k < N; k++) {
                dx      = b_type_d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy      = b_type_d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz      = b_type_d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                pot_vec.store_unaligned(&(ptr[j + k * M]));
            }
        }

        for (std::size_t j = vec_size; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->nr] - this->points_source[cols[k] + 0 * this->nc];
                ddy            = this->points_target[rows[j] + 1 * this->nr] - this->points_source[cols[k] + 1 * this->nc];
                ddz            = this->points_target[rows[j] + 2 * this->nr] - this->points_source[cols[k] + 2 * this->nc];
                ptr[j + k * M] = 1 / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    } else if (N >= simd_size) {
        // std::vector<double> tmp(simd_size);
        std::size_t vec_size = N - N % simd_size;
        for (std::size_t k = 0; k < vec_size; k += simd_size) {
            // load
            Xt_vec.load_unaligned(this->points_source.data() + cols[0] + 0 * this->nc + k);
            Yt_vec.load_unaligned(this->points_source.data() + cols[0] + 1 * this->nc + k);
            Zt_vec.load_unaligned(this->points_source.data() + cols[0] + 2 * this->nc + k);

            for (std::size_t j = 0; j < M; j++) {
                dx      = b_type_d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy      = b_type_d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz      = b_type_d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                // pot_vec.store_aligned(tmp.data());
                for (int i = 0; i < simd_size; i++) {
                    ptr[j + (k + i) * M] = pot_vec[i];
                }
            }
        }
        for (std::size_t k = vec_size; k < N; k++) {
            for (std::size_t j = 0; j < M; ++j) {
                ddx            = this->points_target[rows[j] + 0 * this->nr] - this->points_source[cols[k] + 0 * this->nc];
                ddy            = this->points_target[rows[j] + 1 * this->nr] - this->points_source[cols[k] + 1 * this->nc];
                ddz            = this->points_target[rows[j] + 2 * this->nr] - this->points_source[cols[k] + 2 * this->nc];
                ptr[j + k * M] = 1. / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    } else {
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->nr] - this->points_source[cols[k] + 0 * this->nc];
                ddy            = this->points_target[rows[j] + 1 * this->nr] - this->points_source[cols[k] + 1 * this->nc];
                ddz            = this->points_target[rows[j] + 2 * this->nr] - this->points_source[cols[k] + 2 * this->nc];
                ptr[j + k * M] = 1. / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    }
}

template <>
void MyGenerator<std::complex<double>, 2>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, std::complex<double> *ptr) const {
    using b_type_d                  = xsimd::simd_type<double>;
    using b_type_c_d                = xsimd::simd_type<std::complex<double>>;
    constexpr std::size_t simd_size = b_type_d::size;

    b_type_d charge(1. / (4 * 3.141592653589793238463));
    b_type_d cos, sin, r, aux;
    b_type_d Xt_vec, Yt_vec, Zt_vec;
    b_type_c_d pot_vec;
    b_type_d dx, dy, dz;
    double ddx, ddy, ddz;

    double d;
    if (M >= simd_size) {

        std::size_t vec_size = M - M % simd_size;
        for (std::size_t j = 0; j < vec_size; j += simd_size) {

            // load
            Xt_vec.load_unaligned(this->points_target.data() + rows[0] + 0 * this->nr + j);
            Yt_vec.load_unaligned(this->points_target.data() + rows[0] + 1 * this->nr + j);
            Zt_vec.load_unaligned(this->points_target.data() + rows[0] + 2 * this->nr + j);

            for (std::size_t k = 0; k < N; k++) {
                dx  = b_type_d(this->points_source[cols[0] + k + this->nc * 0]) - Xt_vec;
                dy  = b_type_d(this->points_source[cols[0] + k + this->nc * 1]) - Yt_vec;
                dz  = b_type_d(this->points_source[cols[0] + k + this->nc * 2]) - Zt_vec;
                r   = sqrt(dx * dx + dy * dy + dz * dz);
                aux = charge / r;
                xsimd::sincos(freq * r, sin, cos);
                pot_vec = b_type_c_d(aux * cos, aux * sin);
                pot_vec.store_unaligned(&(ptr[j + k * M]));
                // for (int i = 0; i < simd_size; i++) {
                //     ptr[j + i + k * M] = pot_vec[i];
                // }
            }
        }
        for (std::size_t j = vec_size; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->nr] - this->points_source[cols[k] + 0 * this->nc];
                ddy            = this->points_target[rows[j] + 1 * this->nr] - this->points_source[cols[k] + 1 * this->nc];
                ddz            = this->points_target[rows[j] + 2 * this->nr] - this->points_source[cols[k] + 2 * this->nc];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }
    } else if (N >= simd_size) {
        std::size_t vec_size = N - N % simd_size;
        for (std::size_t k = 0; k < vec_size; k += simd_size) {

            // load
            Xt_vec.load_unaligned(this->points_source.data() + cols[0] + 0 * this->nc + k);
            Yt_vec.load_unaligned(this->points_source.data() + cols[0] + 1 * this->nc + k);
            Zt_vec.load_unaligned(this->points_source.data() + cols[0] + 2 * this->nc + k);

            for (std::size_t j = 0; j < M; j++) {
                dx  = b_type_d(this->points_target[rows[0] + j + this->nr * 0]) - Xt_vec;
                dy  = b_type_d(this->points_target[rows[0] + j + this->nr * 1]) - Yt_vec;
                dz  = b_type_d(this->points_target[rows[0] + j + this->nr * 2]) - Zt_vec;
                r   = sqrt(dx * dx + dy * dy + dz * dz);
                aux = charge / r;
                xsimd::sincos(freq * r, sin, cos);
                pot_vec = b_type_c_d(aux * cos, aux * sin);

                for (int i = 0; i < simd_size; i++) {
                    ptr[j + (k + i) * M] = pot_vec[i];
                }
            }
        }
        for (std::size_t k = vec_size; k < N; k++) {
            for (std::size_t j = 0; j < M; ++j) {
                ddx            = this->points_target[rows[j] + 0 * this->nr] - this->points_source[cols[k] + 0 * this->nc];
                ddy            = this->points_target[rows[j] + 1 * this->nr] - this->points_source[cols[k] + 1 * this->nc];
                ddz            = this->points_target[rows[j] + 2 * this->nr] - this->points_source[cols[k] + 2 * this->nc];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }
    } else {
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->nr] - this->points_source[cols[k] + 0 * this->nc];
                ddy            = this->points_target[rows[j] + 1 * this->nr] - this->points_source[cols[k] + 1 * this->nc];
                ddz            = this->points_target[rows[j] + 2 * this->nr] - this->points_source[cols[k] + 2 * this->nc];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }
    }
}
// xsimd library (aligned)

template <>
void MyGenerator<double, 3>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const {
    using b_type_d                  = xsimd::simd_type<double>;
    constexpr std::size_t simd_size = b_type_d::size;

    b_type_d charge(1. / (4 * 3.141592653589793238463)), r2;
    b_type_d Xt_vec, Yt_vec, Zt_vec;
    b_type_d pot_vec;
    b_type_d dx, dy, dz;
    double ddx, ddy, ddz;
    if (M >= simd_size) {
        std::size_t vec_start = simd_size - rows[0] % simd_size;
        for (std::size_t j = 0; j < vec_start; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                ptr[j + k * M] = 1 / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }

        std::size_t vec_size = (M - vec_start) - (M - vec_start) % simd_size;
        for (std::size_t j = vec_start; j < vec_size; j += simd_size) {

            // load
            Xt_vec.load_aligned(this->points_target.data() + rows[0] + 0 * this->inct + j);
            Yt_vec.load_aligned(this->points_target.data() + rows[0] + 1 * this->inct + j);
            Zt_vec.load_aligned(this->points_target.data() + rows[0] + 2 * this->inct + j);

            for (std::size_t k = 0; k < N; k++) {
                dx      = b_type_d(this->points_source[cols[0] + k + this->incs * 0]) - Xt_vec;
                dy      = b_type_d(this->points_source[cols[0] + k + this->incs * 1]) - Yt_vec;
                dz      = b_type_d(this->points_source[cols[0] + k + this->incs * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                pot_vec.store_unaligned(&(ptr[j + k * M]));
            }
        }
        for (std::size_t j = vec_size; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                ptr[j + k * M] = 1 / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    } else if (N >= simd_size) {
        // std::vector<double, xsimd::aligned_allocator<double>> tmp(8);
        std::size_t vec_start = simd_size - cols[0] % simd_size;
        for (std::size_t k = 0; k < vec_start; k++) {
            for (std::size_t j = 0; j < M; ++j) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                ptr[j + k * M] = 1. / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }

        std::size_t vec_size = (N - vec_start) - (N - vec_start) % simd_size;
        for (std::size_t k = vec_start; k < vec_size; k += simd_size) {
            // load
            Xt_vec.load_aligned(this->points_source.data() + cols[0] + 0 * this->incs + k);
            Yt_vec.load_aligned(this->points_source.data() + cols[0] + 1 * this->incs + k);
            Zt_vec.load_aligned(this->points_source.data() + cols[0] + 2 * this->incs + k);

            for (std::size_t j = 0; j < M; j++) {
                dx      = b_type_d(this->points_target[rows[0] + j + this->inct * 0]) - Xt_vec;
                dy      = b_type_d(this->points_target[rows[0] + j + this->inct * 1]) - Yt_vec;
                dz      = b_type_d(this->points_target[rows[0] + j + this->inct * 2]) - Zt_vec;
                r2      = dx * dx + dy * dy + dz * dz;
                pot_vec = charge / sqrt(r2);
                // pot_vec.store_aligned(tmp.data());
                for (int i = 0; i < simd_size; i++) {
                    ptr[j + (k + i) * M] = pot_vec[i];
                }
            }
        }
        for (std::size_t k = vec_size; k < N; k++) {
            for (std::size_t j = 0; j < M; ++j) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                ptr[j + k * M] = 1. / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    } else {
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                ptr[j + k * M] = 1. / ((4 * 3.141592653589793238463) * (sqrt(ddx * ddx + ddy * ddy + ddz * ddz)));
            }
        }
    }
}

template <>
void MyGenerator<std::complex<double>, 3>::copy_submatrix(int M, int N, const int *const rows, const int *const cols, std::complex<double> *ptr) const {
    using b_type_d                  = xsimd::simd_type<double>;
    using b_type_c_d                = xsimd::simd_type<std::complex<double>>;
    constexpr std::size_t simd_size = b_type_d::size;

    b_type_d charge(1. / (4 * 3.141592653589793238463));
    b_type_d cos, sin, r, aux;
    b_type_d Xt_vec, Yt_vec, Zt_vec;
    b_type_c_d pot_vec;
    b_type_d dx, dy, dz;
    double d, ddx, ddy, ddz;

    if (M >= simd_size) {

        std::size_t vec_start = simd_size - rows[0] % simd_size;
        for (std::size_t j = 0; j < vec_start; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }

        std::size_t vec_size = (M - vec_start) - (M - vec_start) % simd_size;
        for (std::size_t j = vec_start; j < vec_size; j += simd_size) {

            // load
            Xt_vec.load_aligned(this->points_target.data() + rows[0] + 0 * this->inct + j);
            Yt_vec.load_aligned(this->points_target.data() + rows[0] + 1 * this->inct + j);
            Zt_vec.load_aligned(this->points_target.data() + rows[0] + 2 * this->inct + j);

            for (std::size_t k = 0; k < N; k++) {
                dx  = b_type_d(this->points_source[cols[0] + k + this->incs * 0]) - Xt_vec;
                dy  = b_type_d(this->points_source[cols[0] + k + this->incs * 1]) - Yt_vec;
                dz  = b_type_d(this->points_source[cols[0] + k + this->incs * 2]) - Zt_vec;
                r   = sqrt(dx * dx + dy * dy + dz * dz);
                aux = charge / r;
                xsimd::sincos(freq * r, sin, cos);
                pot_vec = b_type_c_d(aux * cos, aux * sin);
                pot_vec.store_unaligned(&(ptr[j + k * M]));
            }
        }
        for (std::size_t j = vec_size; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }
    } else if (N >= simd_size) {
        std::size_t vec_start = simd_size - cols[0] % simd_size;
        for (std::size_t k = 0; k < vec_start; k++) {
            for (std::size_t j = 0; j < M; ++j) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }

        std::size_t vec_size = (N - vec_start) - (N - vec_start) % simd_size;
        for (std::size_t k = vec_start; k < vec_size; k += simd_size) {
            // load
            Xt_vec.load_aligned(this->points_source.data() + cols[0] + 0 * this->incs + k);
            Yt_vec.load_aligned(this->points_source.data() + cols[0] + 1 * this->incs + k);
            Zt_vec.load_aligned(this->points_source.data() + cols[0] + 2 * this->incs + k);

            for (std::size_t j = 0; j < M; j++) {
                dx  = b_type_d(this->points_target[rows[0] + j + this->inct * 0]) - Xt_vec;
                dy  = b_type_d(this->points_target[rows[0] + j + this->inct * 1]) - Yt_vec;
                dz  = b_type_d(this->points_target[rows[0] + j + this->inct * 2]) - Zt_vec;
                r   = sqrt(dx * dx + dy * dy + dz * dz);
                aux = charge / r;
                xsimd::sincos(freq * r, sin, cos);
                pot_vec = b_type_c_d(aux * cos, aux * sin);
                for (int i = 0; i < simd_size; i++) {
                    ptr[j + (k + i) * M] = pot_vec[i];
                }
            }
        }
        for (std::size_t k = vec_size; k < N; k++) {
            for (std::size_t j = 0; j < M; ++j) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }
    } else {
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N; k++) {
                ddx            = this->points_target[rows[j] + 0 * this->inct] - this->points_source[cols[k] + 0 * this->incs];
                ddy            = this->points_target[rows[j] + 1 * this->inct] - this->points_source[cols[k] + 1 * this->incs];
                ddz            = this->points_target[rows[j] + 2 * this->inct] - this->points_source[cols[k] + 2 * this->incs];
                d              = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                ptr[j + k * M] = (1. / (4 * 3.141592653589793238463)) * exp(std::complex<double>(0, 1) * freq * d) / (d);
            }
        }
    }
}
#endif